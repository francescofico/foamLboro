#%%
import os
import numpy as np
import mmap
import re
import gzip
from contextlib import contextmanager


_QUOTED = r'"(?:\\.|[^"\\])*"'
_COMMENTS = re.compile(_QUOTED + r'|//[^\n]*|/\*.*?\*/', re.DOTALL)
_DELIMITERS = re.compile(_QUOTED + r'|[{}()]')
_KEY = re.compile(r'\s*(' + _QUOTED + r'|[^\s{}();"]+)\s*')
_VALUE_END = re.compile(_QUOTED + r'|[;({]')
_HEADER = re.compile(r'\s*FoamFile\s*\{')
_COMPONENTS = {'scalar': 1, 'vector': 3, 'sphericalTensor': 1,
               'symmTensor': 6, 'tensor': 9}
_NUMBER = r'[+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|inf(?:inity)?|nan)'


@contextmanager
def _mapped_data(path):
    """Map ordinary files; compressed ASCII files are decompressed in memory."""
    with open(path, 'rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
        f.seek(0)
        if compressed:
            with gzip.GzipFile(fileobj=f) as stream:
                yield stream.read()
        elif os.fstat(f.fileno()).st_size:
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as data:
                yield data
        else:
            raise ValueError(f'Empty OpenFOAM file: {path}')


def _block_end(text, start):
    opening = text[start]
    closing = {'{': '}', '(': ')'}[opening]
    depth = 1
    for match in _DELIMITERS.finditer(text, start + 1):
        token = match.group()
        if token == opening:
            depth += 1
        elif token == closing:
            depth -= 1
            if not depth:
                return match.end()
    raise ValueError(f'Unclosed {opening!r} in OpenFOAM data')


def _entries(text):
    """Yield direct dictionary entries, skipping balanced nested values."""
    pos = 0
    while pos < len(text):
        if text[pos].isspace() or text[pos] == ';':
            pos += 1
            continue
        match = _KEY.match(text, pos)
        if not match:
            raise ValueError('Invalid OpenFOAM dictionary entry')
        key = match.group(1).strip('"')
        start = pos = match.end()
        if pos < len(text) and text[pos] == '{':
            pos = _block_end(text, pos)
            yield key, text[start:pos]
            continue
        while True:
            end = _VALUE_END.search(text, pos)
            if not end:
                raise ValueError(f'Missing semicolon for entry {key!r}')
            token = end.group()
            if token == ';':
                yield key, text[start:end.start()].strip()
                pos = end.end()
                break
            pos = _block_end(text, end.start()) if token in ('(', '{') else end.end()


def _entry(text, key, required=True):
    for name, value in _entries(text):
        if name == key:
            return value
    if required:
        raise ValueError(f'OpenFOAM entry {key!r} not found')
    return None


def _ascii_text(data):
    try:
        text = data.decode('utf-8')
    except UnicodeDecodeError as exc:
        raise ValueError('Only ASCII OpenFOAM data are supported; convert binary fields to ASCII') from exc
    text = _COMMENTS.sub(lambda m: m.group() if m.group().startswith('"') else ' ', text)
    header = _HEADER.match(text)
    if header:
        body = text[header.end():_block_end(text, header.end() - 1) - 1]
        if (_entry(body, 'class', False) or '').strip('"') == 'decomposedBlockData':
            raise ValueError('Use the collated reconstruction functions for decomposedBlockData')
        if (_entry(body, 'format', False) or 'ascii').strip('"') != 'ascii':
            raise ValueError('Only ASCII OpenFOAM data are supported; convert binary fields to ASCII')
    if '\x00' in text:
        raise ValueError('Binary data found; only ASCII OpenFOAM data are supported')
    return text


def _read_text(path):
    with _mapped_data(path) as data:
        return _ascii_text(data[:])


def _without_header(text):
    header = _HEADER.match(text)
    return text[_block_end(text, header.end() - 1):].strip() if header else text.strip()


def _counted_list(text):
    match = re.match(r'\s*(?:List<\w+>\s+)?(\d+)\s*', text)
    if not match:
        raise ValueError('Missing OpenFOAM list count')
    count, start = int(match.group(1)), match.end()
    if count == 0 and text[start:].strip() in ('', ';'):
        return 0, '', False
    if start == len(text) or text[start] not in '({':
        raise ValueError('Missing delimiters for OpenFOAM list')
    end = _block_end(text, start)
    if text[end:].strip() not in ('', ';'):
        raise ValueError('Unexpected data after OpenFOAM list')
    return count, text[start + 1:end - 1], text[start] == '{'


def _numeric_values(text, components=1, tuples=False, integer=False):
    number = r'[+-]?\d+' if integer else _NUMBER
    element = rf'\(\s*{number}(?:\s+{number}){{{components - 1}}}\s*\)' if tuples else number
    pattern = re.compile(element, re.IGNORECASE)
    pos = count = 0
    for match in pattern.finditer(text):
        if text[pos:match.start()].strip():
            raise ValueError('Invalid ASCII numeric data or component count')
        pos, count = match.end(), count + 1
    if text[pos:].strip():
        raise ValueError('Invalid ASCII numeric data or component count')
    if not count:
        return np.empty((0, components), dtype=np.int64 if integer else np.float64)
    if integer:
        arr = np.fromiter((int(m.group()) for m in pattern.finditer(text)), dtype=np.int64, count=count)
    else:
        arr = np.fromstring(text.translate(str.maketrans('()', '  ')), sep=' ', dtype=np.float64)
    if arr.size != count * components:
        raise ValueError('Invalid ASCII numeric data or component count')
    return arr.reshape(count, components)


def _list_values(text, components=1, tuples=False, integer=False, expected_count=None):
    count, body, repeated = _counted_list(text)
    if expected_count is not None and count != expected_count:
        raise ValueError(f'Expected {expected_count} entries, but the file declares {count}')
    values = _numeric_values(body, components, tuples, integer)
    if len(values) != (1 if repeated else count):
        raise ValueError(f'List count mismatch: declared {count}, read {len(values)} entries')
    return np.repeat(values, count, axis=0) if repeated else values


def _field_values(text, key, kind=None, count=None):
    value = _entry(text, key)
    nonuniform = re.match(r'nonuniform\s+List<(\w+)>\s*', value)
    if nonuniform:
        actual_kind = nonuniform.group(1)
        if actual_kind not in _COMPONENTS or (kind and kind != actual_kind):
            raise ValueError(f'Unexpected field type {actual_kind!r} for {key}')
        return _list_values(value[nonuniform.end():], _COMPONENTS[actual_kind],
                            actual_kind != 'scalar', expected_count=count)
    uniform = re.match(r'uniform\s+', value)
    if not uniform:
        raise ValueError(f'Expected uniform or nonuniform data for {key}')
    body = value[uniform.end():].strip()
    tuples = body.startswith('(')
    components = _COMPONENTS[kind] if kind else 1
    if kind is None and tuples:
        end = _block_end(body, 0)
        components = len(_numeric_values(body[1:end - 1]))
        if components not in (1, 3, 6, 9):
            raise ValueError('Unsupported uniform field component count')
    values = _numeric_values(body, components, tuples=(kind != 'scalar' if kind else tuples))
    if len(values) != 1:
        raise ValueError(f'Expected one uniform value for {key}')
    return np.repeat(values, 1 if count is None else count, axis=0)


def _processor_blocks(path):
    """Decode one ASCII processor payload at a time using its byte count."""
    marker = re.compile(rb'(?m)^[ \t]*(?://[ \t]*)?processor(\d+)\s+(?:List<char>\s+)?(\d+)\s*', re.IGNORECASE)
    with _mapped_data(path) as data:
        pos, seen = 0, set()
        while True:
            match = marker.search(data, pos)
            if not match:
                break
            if not seen:
                # The outer container is binary even for ASCII payloads.
                # v2606 records the payload format separately as data.format.
                prefix = data[:match.start()].decode('utf-8')
                prefix = _COMMENTS.sub(lambda m: m.group() if m.group().startswith('"') else ' ', prefix)
                header = _HEADER.match(prefix)
                if header:
                    body = prefix[header.end():_block_end(prefix, header.end() - 1) - 1]
                    fmt = _entry(body, 'data.format', False)
                    if fmt and fmt.strip('"') != 'ascii':
                        raise ValueError('Only ASCII OpenFOAM payloads are supported; convert binary fields to ASCII')
            processor, size = int(match.group(1)), int(match.group(2))
            if processor in seen:
                raise ValueError(f'Duplicate processor block {processor}')
            seen.add(processor)
            start = match.end()
            if size:
                end = start + 1 + size
                if data[start:start + 1] != b'(' or data[end:end + 1] != b')':
                    raise ValueError(f'Truncated or invalid byte count for processor {processor}')
                yield processor, _ascii_text(data[start + 1:end])
                pos = end + 1
            else:
                if data[start:start + 2] == b'()':
                    start += 2
                yield processor, ''
                pos = start
        if not seen:
            raise ValueError('No processor blocks found; expected collated OpenFOAM data')


def _validate_addressing(proc_global_cells):
    if not proc_global_cells:
        raise ValueError('No processor addressing supplied')
    mapping = {}
    for processor, cells in proc_global_cells.items():
        indices = np.asarray(cells)
        if not isinstance(processor, (int, np.integer)) or processor < 0:
            raise ValueError('Processor IDs must be non-negative integers')
        if indices.ndim != 1 or (indices.size and indices.dtype.kind not in 'iu'):
            raise ValueError(f'Processor {processor} addressing must be a one-dimensional integer array')
        mapping[processor] = indices
    ncells = sum(indices.size for indices in mapping.values())
    covered = np.zeros(ncells, dtype=bool)
    for processor, indices in mapping.items():
        if np.any(indices < 0) or np.any(indices >= ncells):
            raise ValueError('Cell addressing must cover global indices 0 through Ncells-1')
        indices = mapping[processor] = indices.astype(np.int64, copy=False)
        if np.unique(indices).size != indices.size or np.any(covered[indices]):
            raise ValueError('Duplicate global cell indices in processor addressing')
        covered[indices] = True
    return mapping, ncells


def _check_processors(seen, mapping):
    if seen != set(mapping):
        raise ValueError(f'Missing processor blocks: {sorted(set(mapping) - seen)}')


def _boundary_patches(path):
    count, body, repeated = _counted_list(_without_header(_read_text(path)))
    patches = list(_entries(body))
    if repeated or len(patches) != count or len({name for name, _ in patches}) != count:
        raise ValueError('Boundary patch count mismatch or duplicate patch names')
    for name, value in patches:
        if not value.startswith('{'):
            raise ValueError(f'Expected dictionary for patch {name!r}')
        yield name, value[1:-1]


def find_latesttime(path):
    dir_list = os.listdir(path)
    time_list = []

    for directory in dir_list:
        try:
            time = float(directory)
        except ValueError:
            continue
        if np.isfinite(time) and os.path.isdir(os.path.join(path, directory)):
            time_list.append(directory)
    if not time_list:
        raise ValueError(f'No numeric time directories found in {path}')
    time_list.sort(key=float)
    return time_list[-1]

def read_field_mmap(file_path):
    """Read an ASCII internalField: scalars (N,), vectors/tensors (N, components).

    A uniform field returns one value/row for broadcasting: the field file
    alone does not specify its cell count. ASCII .gz files are also accepted.
    """
    data = _field_values(_read_text(file_path), 'internalField')
    return data[:, 0] if data.shape[1] == 1 else data


def read_boundary_field_vector(filepath, patch_name):
    """
    Reads the value entry of an ASCII vector boundary field.

    Args:
        filepath (str): Path to the OpenFOAM field file (e.g., 0/wallShearStress).
        patch_name (str): Name of the boundary patch (e.g., 'bottom').

    Returns:
        np.ndarray: Array of shape (N, 3). Uniform values return one row
        for broadcasting; the field file does not specify the patch size.
    """
    boundary = _entry(_read_text(filepath), 'boundaryField')
    if not boundary.startswith('{'):
        raise ValueError('Expected a boundaryField dictionary')
    patch_block = _entry(boundary[1:-1], patch_name)
    if not patch_block.startswith('{'):
        raise ValueError(f'Expected a dictionary for patch {patch_name!r}')
    return _field_values(patch_block[1:-1], 'value', 'vector')

def parse_cellProcAddressing(path):
    """
    Parse collated ASCII cellProcAddressing to get global indices per processor.
    Returns a dict: proc_global_cells[proc_id] = array of global indices
    """
    proc_global_cells = {}
    for proc_id, block in _processor_blocks(path):
        indices = _list_values(_without_header(block) or '0()', integer=True)[:, 0]
        if np.any(indices < 0):
            raise ValueError(f'Negative global cell index in processor {proc_id}')
        proc_global_cells[proc_id] = indices

    return proc_global_cells


def reconstruct_vector_mmap(collated_path, proc_global_cells):
    """
    Reconstruct a collated ASCII vector field in global OpenFOAM cell order.
    Uniform processor fields are expanded using their addressing counts.
    ---
    collated only
    """
    proc_global_cells, Ncells = _validate_addressing(proc_global_cells)
    global_U = np.empty((Ncells, 3), dtype=np.float64)
    seen = set()
    for proc_id, block_text in _processor_blocks(collated_path):
        if proc_id not in proc_global_cells:
            raise ValueError(f'No addressing for processor {proc_id}')
        seen.add(proc_id)
        indices = proc_global_cells[proc_id]
        vectors = (_field_values(block_text, 'internalField', 'vector', len(indices))
                   if block_text.strip() else np.empty((0, 3)))
        if len(vectors) != len(indices):
            raise ValueError(f'Processor {proc_id} local cell count mismatch')
        global_U[indices, :] = vectors
    _check_processors(seen, proc_global_cells)
    return global_U


def reconstruct_scalar_field(collated_path, proc_global_cells):
    """
    Reconstruct scalar OpenFOAM field from collated decomposed format.
    Returns shape (Ncells, 1); uniform processor fields are expanded.
    ---
    collated only

    """
    proc_global_cells, Ncells = _validate_addressing(proc_global_cells)
    global_field = np.empty((Ncells, 1), dtype=np.float64)
    seen = set()
    for proc_id, block_text in _processor_blocks(collated_path):
        if proc_id not in proc_global_cells:
            raise ValueError(f'No addressing for processor {proc_id}')
        seen.add(proc_id)
        indices = proc_global_cells[proc_id]
        arr = (_field_values(block_text, 'internalField', 'scalar', len(indices))
               if block_text.strip() else np.empty((0, 1)))
        if len(indices) != arr.shape[0]:
            raise ValueError(
                f"Processor {proc_id}: cell count mismatch {len(indices)} vs {arr.shape[0]}"
            )
        global_field[indices, 0] = arr[:, 0]
    _check_processors(seen, proc_global_cells)
    return global_field


def reconstruct_symmTensor_field(collated_path, proc_global_cells):
    """
    Reconstruct symmetric tensor OpenFOAM field from collated decomposed format.
    Symmetric tensors have 6 components: (xx, xy, xz, yy, yz, zz).
    Uniform processor fields are expanded using their addressing counts.
    
    Parameters
    ----------
    collated_path : str
        Path to collated decomposed field (ASCII format).
    proc_global_cells : dict
        Mapping proc_id -> global cell indices (from parse_cellProcAddressing).
    
    Returns
    -------
    global_field : np.ndarray
        Array of shape (Ncells, 6)
    """
    nComp = 6  # symmetric tensor
    proc_global_cells, Ncells = _validate_addressing(proc_global_cells)
    global_field = np.empty((Ncells, nComp), dtype=np.float64)
    seen = set()
    for proc_id, block_text in _processor_blocks(collated_path):
        if proc_id not in proc_global_cells:
            raise ValueError(f'No addressing for processor {proc_id}')
        seen.add(proc_id)
        indices = proc_global_cells[proc_id]
        arr = (_field_values(block_text, 'internalField', 'symmTensor', len(indices))
               if block_text.strip() else np.empty((0, nComp)))
        if len(indices) != arr.shape[0]:
            raise ValueError(
                f"Processor {proc_id}: cell count mismatch {len(indices)} vs {arr.shape[0]}"
            )
        global_field[indices, :] = arr
    _check_processors(seen, proc_global_cells)
    return global_field



def list_wall_patches(boundary_file, types=('wall', 'mappedWall')):
    """
    List patches of specified types (e.g., wall, mappedWall) from constant/polyMesh/boundary.

    Parameters:
        boundary_file (str): Path to constant/polyMesh/boundary
        types (tuple): Patch types to include

    Returns:
        List[str]: Names of patches with specified types
    """
    patches = []
    for patch_name, patch in _boundary_patches(boundary_file):
        patch_type = _entry(patch, 'type').strip('"')
        if patch_type in types:
            patches.append(patch_name)
    return patches

def read_boundary_faces(path):
    """
    read boundary faces from OpenFOAM boundary file
    Args:
        path (str): Path to the OpenFOAM boundary file
    Returns:
        dict: Dictionary containing patch names as keys and their nFaces and startFace as values
    """
    patch_info = {}
    for patch_name, patch in _boundary_patches(path):
        nFaces = int(_entry(patch, 'nFaces'))
        startFace = int(_entry(patch, 'startFace'))
        if nFaces < 0 or startFace < 0:
            raise ValueError(f'Negative face count or startFace for patch {patch_name!r}')
        patch_info[patch_name] = {'nFaces': nFaces, 'startFace': startFace}

    return patch_info




def detect_regions(collated_time_path):
    """Return region subdirectories, excluding standard OpenFOAM metadata."""
    dirs = [
        d for d in os.listdir(collated_time_path)
        if d not in ('uniform', 'polyMesh', 'lagrangian')
        and os.path.isdir(os.path.join(collated_time_path, d))
    ]
    return dirs  # may be empty → single-region case



def read_owner_file_mmap(path):
    """Read an ASCII owner labelList, retaining 64-bit cell indices."""
    data = _without_header(_read_text(path))
    owner_array = _list_values(data, integer=True)[:, 0]
    if np.any(owner_array < 0):
        raise ValueError('Negative cell index in owner data')
    return owner_array


