# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import concurrent.futures

def system_compress(file_path):
    if not os.path.exists(file_path):
        return f"[SKIP] File not found: {file_path}"
    if file_path.endswith('.gz'):
        return f"[SKIP] Already compressed: {file_path}"

    use_pigz = shutil.which("pigz") is not None
    tool = "pigz" if use_pigz else "gzip"
    cmd = [tool, "-f", file_path]

    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        return f"[COMPRESS] Success: {file_path} -> {file_path}.gz (Engine: {tool})"
    except subprocess.CalledProcessError as e:
        return f"[ERROR] Compression failed {file_path}: {e}"

def compress_inputs(file_list, logger=None):
    valid_files = [f for f in file_list if f]

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(valid_files)) as executor:
        results = executor.map(system_compress, valid_files)
        for res in results:
            if logger:
                logger.info(f"         {res}")