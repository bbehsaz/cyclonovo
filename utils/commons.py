import xml.etree.ElementTree
from distutils import dir_util
from os.path import join, isfile, isdir, basename


def parse_params_xml(fpath):
    import xml.etree.ElementTree
    from distutils import dir_util
    from os.path import join, isfile, isdir, basename


    with open(fpath) as f:
        content = f.read()
    params = dict()
    file_mapping = dict()
    for e in xml.etree.ElementTree.fromstring(content).findall('parameter'):
        if e.attrib['name'] == 'upload_file_mapping':
            fname, real_fname = e.text.split('|')[0:2]
            file_mapping[basename(fname)] = real_fname
            continue
        if e.text == 'on':
            value = True
        elif e.text == 'off':
            value = False
        else:
            value = e.text
        params[e.attrib['name']] = value
    return params, file_mapping
    
def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.
    Parameters
    """
    from distutils import dir_util
    from os.path import join, isfile, isdir, basename
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg




def get_original_fpath(fpath, file_mapping=None):
    import xml.etree.ElementTree
    from distutils import dir_util
    from os.path import join, isfile, isdir, basename
    return file_mapping[basename(fpath)] if file_mapping and basename(fpath) in file_mapping else fpath

