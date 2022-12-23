def get_free_memory():
    """ Return available memory in kB"""
    import re
    with open("/proc/meminfo") as f:
        meminfo = f.read()
        matched = re.search(r"MemFree:\s+(\d+)", meminfo)

    return int(matched.groups()[0])


available_mem = get_free_memory()

