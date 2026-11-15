"""Regression checks for the reader's public functions """
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


spec = importlib.util.spec_from_file_location('reader', Path(__file__).resolve().parents[1]/'src/reader.py')
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)


def header(kind='volScalarField', fmt='ascii'):
    return f'FoamFile {{ version 2.0; format {fmt}; class {kind}; object test; }}\n'


def collated(blocks, prefix='processor'):
    data = header('decomposedBlockData', 'binary').encode()
    for processor, payload in blocks:
        payload = payload.encode() if isinstance(payload, str) else payload
        data += f'// {prefix}{processor}\n{len(payload)}\n'.encode()
        data += b'('+payload+b')\n' if payload else b'\n'
    return data


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='foam-reader-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def file(self, content, name='field'):
        path = self.root/name
        path.write_bytes(content.encode() if isinstance(content, str) else content)
        return path

    def test_internal_fields_and_stable_shapes(self):
        for body, expected in (
                ('nonuniform List<scalar> 3(1 -2.5 3e-4)', [1, -2.5, 3e-4]),
                ('nonuniform List<scalar>\n1\n(\n7\n)', [7]),
                ('nonuniform List<scalar> 0()', []),
                ('nonuniform List<scalar> 0( \n\t )', []),
                ('nonuniform List<scalar> 3{7}', [7, 7, 7]),
                ('nonuniform List<vector> 2((1 2 3)(4 5 6))', [[1, 2, 3], [4, 5, 6]]),
                ('nonuniform List<vector>\n1\n(\n(1 2 3)\n)', [[1, 2, 3]]),
                ('nonuniform List<vector> 0()', np.empty((0, 3))),
                ('nonuniform List<vector> 0(/* empty */\n)', np.empty((0, 3))),
                ('nonuniform List<symmTensor> 1((1 2 3 4 5 6))', [[1, 2, 3, 4, 5, 6]]),
                ('uniform (4 5 6)', [[4, 5, 6]])):
            with self.subTest(body=body):
                actual = reader.read_field_mmap(self.file('internalField '+body+';'))
                np.testing.assert_array_equal(actual, expected)
                self.assertEqual(actual.shape, np.asarray(expected).shape)

    def test_uniform_internal_field_does_not_read_boundary_values(self):
        path = self.file(header()+'internalField uniform 9;\nboundaryField { wall { '
                         'value nonuniform List<scalar> 2(100 200); } }')
        np.testing.assert_array_equal(reader.read_field_mmap(path), [9])
        with self.assertRaisesRegex(ValueError, 'internalField'):
            reader.read_field_mmap(self.file('boundaryField { wall { internalField uniform 99; } }'))

    def test_comments_and_nested_exact_patch_lookup(self):
        path = self.file(header('volVectorField')+'// internalField uniform (99 99 99);\n'
            'internalField uniform (0 0 0); boundaryField { '
            'frontwall { value uniform (1 2 3); } '
            'wall { options { value uniform (99 99 99); note "} // quoted"; } '
            '/* value uniform (8 8 8); } */ value nonuniform List<vector> 1((7 8 9)); } '
            'wallx1 { value uniform (10 11 12); } "wall.1" { value uniform (13 14 15); } }')
        np.testing.assert_array_equal(reader.read_field_mmap(path), [[0, 0, 0]])
        np.testing.assert_array_equal(reader.read_boundary_field_vector(path, 'wall'), [[7, 8, 9]])
        np.testing.assert_array_equal(reader.read_boundary_field_vector(path, 'wall.1'), [[13, 14, 15]])
        with self.assertRaisesRegex(ValueError, 'not found'):
            reader.read_boundary_field_vector(path, 'missing')

    def test_malformed_fields_fail(self):
        for body in ('nonuniform List<scalar> 3(1 2)',
                     'nonuniform List<scalar> 2(1 2 junk)',
                     'nonuniform List<vector> 2((1 2)(3 4 5 6))',
                     'nonuniform List<vector> 1((1 2 3)',
                     'uniform (1 2)', 'uniform 1 2'):
            with self.subTest(body=body), self.assertRaises(ValueError):
                reader.read_field_mmap(self.file('internalField '+body+';'))

    def test_collated_addressing_and_all_field_types(self):
        mapping_data = [(1, '2(3 1)'), (0, header('labelList')+'2(2 0)'), (2, '0()')]
        for prefix in ('processor', 'Processor'):
            mapping = reader.parse_cellProcAddressing(self.file(collated(mapping_data, prefix)))
            for pid, indices in ((0, [2, 0]), (1, [3, 1]), (2, [])):
                np.testing.assert_array_equal(mapping[pid], indices)
            for function, kind, first, second, expected in (
                    (reader.reconstruct_scalar_field, 'scalar', 'uniform 10',
                     'nonuniform List<scalar> 2(30 40)', [[10], [40], [10], [30]]),
                    (reader.reconstruct_vector_mmap, 'vector',
                     'nonuniform List<vector> 2((21 22 23) (1 2 3))', 'uniform (7 8 9)',
                     [[1, 2, 3], [7, 8, 9], [21, 22, 23], [7, 8, 9]]),
                    (reader.reconstruct_symmTensor_field, 'symmTensor',
                     'nonuniform List<symmTensor> 2((1 2 3 4 5 6)(7 8 9 10 11 12))',
                     'uniform (13 14 15 16 17 18)',
                     [[7, 8, 9, 10, 11, 12], [13, 14, 15, 16, 17, 18],
                      [1, 2, 3, 4, 5, 6], [13, 14, 15, 16, 17, 18]])):
                with self.subTest(prefix=prefix, field=kind):
                    blocks = [(1, 'internalField '+second+';'),
                              (0, header()+'internalField '+first+';'),
                              (2, f'internalField nonuniform List<{kind}> 0();')]
                    path = self.file(collated(blocks, prefix))
                    np.testing.assert_array_equal(function(path, mapping), expected)
                    self.assertEqual(function(path, mapping).shape, np.asarray(expected).shape)

    def test_reconstruction_rejects_incomplete_or_duplicate_processors(self):
        for function, value in ((reader.reconstruct_scalar_field, '1'),
                                (reader.reconstruct_vector_mmap, '(1 2 3)'),
                                (reader.reconstruct_symmTensor_field, '(1 2 3 4 5 6)')):
            payload = 'internalField uniform '+value+';'
            for blocks, error in (([(0, payload)], 'Missing processor'),
                                  ([(0, payload), (0, payload)], 'Duplicate processor'),
                                  ([(0, payload), (2, payload)], 'No addressing')):
                with self.subTest(function=function.__name__, blocks=blocks), self.assertRaisesRegex(ValueError, error):
                    function(self.file(collated(blocks)), {0: [0], 1: [1]})

    def test_invalid_addressing_is_rejected(self):
        path = self.file(collated([(0, 'internalField uniform 1;')]))
        for mapping in ({}, {0: [0, 0]}, {0: [0], 1: [0]}, {0: [-1, 0]},
                        {0: [0, 2]}, {0: [0.5]}, {0: [[0]]}, {-1: [0]}):
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                reader.reconstruct_scalar_field(path, mapping)
        with self.assertRaisesRegex(ValueError, 'count mismatch'):
            reader.parse_cellProcAddressing(self.file(collated([(0, '3(0 1)')])))
        with self.assertRaisesRegex(ValueError, 'declares 1'):
            reader.reconstruct_scalar_field(
                self.file(collated([(0, 'internalField nonuniform List<scalar> 1(7);')])), {0: [0, 1]})

    def test_zero_cell_processors(self):
        mapping = reader.parse_cellProcAddressing(self.file(collated([(0, '')])))
        self.assertEqual(mapping[0].shape, (0,))
        for function, shape in ((reader.reconstruct_scalar_field, (0, 1)),
                                (reader.reconstruct_vector_mmap, (0, 3)),
                                (reader.reconstruct_symmTensor_field, (0, 6))):
            self.assertEqual(function(self.file(collated([(0, '')])), mapping).shape, shape)

    def test_collated_byte_counts_and_embedded_markers(self):
        payload = '/*\n// processor99\n1(x)\n*/ internalField uniform 7;'
        data = collated([(0, payload)])
        np.testing.assert_array_equal(reader.reconstruct_scalar_field(self.file(data), {0: [0]}), [[7]])
        for invalid in (data[:-2], data.replace(f'\n{len(payload)}\n'.encode(), b'\n9999\n', 1)):
            with self.assertRaisesRegex(ValueError, 'byte count'):
                reader.reconstruct_scalar_field(self.file(invalid), {0: [0]})
        with self.assertRaisesRegex(ValueError, 'No processor blocks'):
            reader.parse_cellProcAddressing(self.file(header('labelList')+'2(0 1)'))

    def test_binary_fields_are_rejected(self):
        for payload in (header(fmt='binary')+'internalField uniform 3;',
                        header(fmt='"binary"')+'internalField uniform 3;',
                        b'internalField nonuniform List<scalar> 1(\x00\x00\xff);'):
            with self.subTest(payload=payload), self.assertRaisesRegex(ValueError, 'ASCII'):
                reader.read_field_mmap(self.file(payload))
            with self.subTest(collated=payload), self.assertRaisesRegex(ValueError, 'ASCII'):
                reader.reconstruct_scalar_field(self.file(collated([(0, payload)])), {0: [0]})

    def test_collated_payload_format_is_distinct_from_container_format(self):
        data = collated([(0, 'internalField uniform 3;')])
        for fmt in ('ascii', 'binary'):
            path = self.file(data.replace(b'object test;', f'object test; data.format {fmt};'.encode(), 1))
            if fmt == 'ascii':
                np.testing.assert_array_equal(reader.reconstruct_scalar_field(path, {0: [0]}), [[3]])
            else:
                with self.assertRaisesRegex(ValueError, 'ASCII'):
                    reader.reconstruct_scalar_field(path, {0: [0]})

    def test_compressed_ascii_files(self):
        path = self.file(gzip.compress(b'internalField nonuniform List<scalar> 2(4 5);'), 'T.gz')
        np.testing.assert_array_equal(reader.read_field_mmap(path), [4, 5])
        path = self.file(gzip.compress(collated([(0, 'internalField uniform 6;')])), 'collated.gz')
        np.testing.assert_array_equal(reader.reconstruct_scalar_field(path, {0: [0, 1]}), [[6], [6]])

    def test_boundary_mesh_helpers(self):
        path = self.file(header('polyBoundaryMesh')+'3(\n'
            'fluid-wall { startFace 20; options { type patch; } nFaces 10; type mappedWall; }\n'
            '"wall.1"\n{\ntype wall; inGroups 1(wall); nFaces 4; startFace 30;\n}\n'
            'outlet { type patch; nFaces 2; startFace 34; }\n)')
        self.assertEqual(reader.list_wall_patches(path), ['fluid-wall', 'wall.1'])
        self.assertEqual(reader.read_boundary_faces(path), {
            'fluid-wall': {'nFaces': 10, 'startFace': 20},
            'wall.1': {'nFaces': 4, 'startFace': 30},
            'outlet': {'nFaces': 2, 'startFace': 34}})
        with self.assertRaisesRegex(ValueError, 'count mismatch'):
            reader.read_boundary_faces(self.file('2(wall { type wall; nFaces 1; startFace 0; })'))

    def test_owner_counts_and_large_indices(self):
        for body, expected in (('3(0 1 3000000000)', [0, 1, 3000000000]),
                               ('3{4}', [4, 4, 4]), ('0()', [])):
            np.testing.assert_array_equal(reader.read_owner_file_mmap(self.file(header('labelList')+body)), expected)
        for body in ('3(0 1)', '2(-1 0)', '2(0 1 bad)'):
            with self.subTest(body=body), self.assertRaises(ValueError):
                reader.read_owner_file_mmap(self.file(header('labelList')+body))

    def test_latest_time_requires_finite_numeric_directories(self):
        for name in ('constant', 'system', '1.5', '1e2', 'nan', 'inf'):
            (self.root/name).mkdir()
        (self.root/'1000').write_text('numeric filename')
        self.assertEqual(reader.find_latesttime(self.root), '1e2')
        with self.assertRaisesRegex(ValueError, 'No numeric time directories'):
            reader.find_latesttime(self.root/'constant')

    def test_region_detection_excludes_metadata(self):
        for name in ('fluid', 'solid', 'uniform', 'polyMesh', 'lagrangian'):
            (self.root/name).mkdir()
        (self.root/'U').write_text('field')
        self.assertEqual(set(reader.detect_regions(self.root)), {'fluid', 'solid'})


if __name__ == '__main__':
    unittest.main()
