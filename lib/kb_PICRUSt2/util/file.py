import re
import os

from .cli import gunzip

def get_numbered_duplicate(names, q):
    if q not in names:
        return q
    qregex = '^' + re.escape(q) + r'( \([1-9]\d*\))?$' # the_attr, the_attr (1), the_attr (2) ...
    nums = []
    for name in names:
        if re.match(qregex, name) and name != q:
            num = int(name.split('(')[-1][:-1])
            nums.append(num)
    i = 1
    while True:
        if i not in nums:
            break
        i = i + 1
    return q + ' (%d)' % i


def gunzip_out(src_flpth_l, dst_dir, gz_ext='.gz'):
    '''
    Given list of gzipped files, gunzip them to another location
    '''
    if not os.path.exists(dst_dir):
        os.mkdir(dst_dir)

    dst_flpth_l = []
    for i, src_flpth in enumerate(src_flpth_l):
        subdir = os.path.join(dst_dir, str(i)) # create subdirectories 0, 1, 2, ... since src files can have same name
        os.mkdir(subdir)
        dst_flpth = os.path.join(subdir, os.path.basename(src_flpth[:-len(gz_ext)]))

        gunzip(
            src_flpth,
            dst_flpth
        )
        
        dst_flpth_l.append(dst_flpth)

    return dst_flpth_l
