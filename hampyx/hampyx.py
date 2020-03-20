import os
import tempfile
import subprocess
import xml.etree.ElementTree as et
import numpy as np


from functools import wraps
def icy(cls):
    """
    A class decorator for freezing attributes
    """
    cls.__frozen = False
    def frozensetattr(self, key, value):
        if self.__frozen and not hasattr(self, key):
            print("Class {} is frozen. Cannot set {} = {}"
                  .format(cls.__name__, key, value))
        else:
            object.__setattr__(self, key, value)
    def init_decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            func(self, *args, **kwargs)
            self.__frozen = True
        return wrapper
    cls.__setattr__ = frozensetattr
    cls.__init__ = init_decorator(cls.__init__)
    return cls


@icy
class Hampyx(object):
    """
    methods:
    
    # Import class
    
    In []: import hampyx as hpx
    
    # Initialize object
    
    In []: object = hpx.Hampyx (xml_path=<xml file path>, exe_path=<executable path>)
    
    # hammurabiX executable path,
    # by default,
    # is searched from the environment varialbe PATH,
    # while xml file path,
    # by default,
    # is searched in the current working directory './'.
    
    # Modify parameter value from base xml file to temp XML file
    
    In []: object.mod_par (keychain=['key1','key2',...], attrib={'tag':'content'})
    
    # Add new parameter
    
    In []: object.add_par (keychain=['key1','key2',...], subkey='keyfinal', attrib={'tag':'content'})
    
    # Delete parameter
    
    In []: object.del_par (keychain=['key1','key2',...])
    
    # If additional argument opt='all', then all matching parameters will be deleted.
    # The strings 'key1', 'key2',
    # etc represent the path to the desired parameter, going through the XML.
    # The "tag" is the label for the parameter: eg. "Value" or "cue" or "type".
    # The "content" is the content under the tag: eg. the string for the tag "filename".
    
    # Look through the parameter tree in python
    
    In []: object.print_par(keychain=['key1','key2',...])
    
    # This will return the current value of the parameter in the XML,
    # associated with the path "key1/key2/.../keyfinal/".
    
    # Run the executable
    
    In []: object(raiseonerr=True/False)
    
    # If additional verbose=True (by default is False),
    # run.log and err.log will be dumped to disk,
    # notice that dumping logs is not thread safe, use quiet mode in threading.
    
    # After this main routine,
    # object.sim_map will be filled with simulation outputs from hammurabiX.
    # The structure of object.sim_map contains arrays under entries:
    # (we give up nested dict structure for the convenience of Bayesian analysis)
    # object.sim_map[('sync',str(freq),str(Nside),'I')] # synchrotron intensity map at 'frequency'
    # object.sim_map[('sync',str(freq),str(Nside),'Q')] # synchrotron Q map at 'frequency'
    # object.sim_map[('sync',str(freq),str(Nside),'U')] # synchrotron U map at 'frequency'
    # object.sim_map[('fd','nan',str(Nside),'nan')] # Faraday depth map
    # object.sim_map[('dm','nan',str(Nside),'nan')] # dispersion measure map
    """
    
    def __init__(self, xml_path=None, exe_path=None):
        """
        Hampyx initialization,
        default executable path is None, we will search users' environment,
        default parameter file path is None, we will search current working directory,
        otherwise, ABSOLUTE paths are required.
        
        Parameters
        ----------
        
        xml_path : string
            XML parameter file path
            
        exe_path : string
            hammurabi C executable path
        """
        # current working directory
        self.wk_dir = None
        # encapsulated below
        self.exe_path = exe_path
        self.xml_path = xml_path
        # assign tmp file name by base file name
        self.temp_file = self._base_file
        # read from base parameter file
        self.tree = et.parse(self._base_file)
        # simulation output
        self.sim_map_name = {}
        self.sim_map = {}
        # switches
        self._do_sync = False
        self._do_dm = False
        self._do_fd = False
    
    @property
    def exe_path(self):
        return self._exe_path
    
    @property
    def xml_path(self):
        return self._xml_path
    
    @exe_path.setter
    def exe_path(self, exe_path):
        """
        By default hammurabiX executable "hamx" should be available in PATH.
        """
        if exe_path is None:  # search sys environ
            self._exe_path = None
            env = os.environ.get('PATH').split(os.pathsep)
            for dir in env:  # get the first match that exists
                if os.path.isfile(os.path.join(dir, 'hamx')):
                    self._exe_path = os.path.join(dir, 'hamx')
                    break
            assert (self._exe_path is not None)
        else:  # if given
            assert isinstance(exe_path, str)
            self._exe_path = os.path.abspath(exe_path)
        assert (os.path.isfile(self._exe_path))
        self._executable = self._exe_path

    @xml_path.setter
    def xml_path(self, xml_path):
        """
        By default an "*.xml" file should be accessible in current working directory,
        to be specific, the return of os.getcwd().
        """
        if xml_path is None:
            self._xml_path = None
            cnddt = [s for s in os.listdir(self._wk_dir) if 'xml' in s]
            for match in cnddt:  # get the first match that exists
                if os.path.isfile(os.path.join(self._wk_dir, match)):
                    self._xml_path = os.path.join(self._wk_dir, match)
                    break
            assert (self._xml_path is not None)
        else:
            assert isinstance(xml_path, str)
            self._xml_path = os.path.abspath(xml_path)
        assert os.path.isfile(self._xml_path)
        self._base_file = self._xml_path

    @property
    def wk_dir(self):
        return self._wk_dir

    @wk_dir.setter
    def wk_dir(self, wk_dir):
        """
        The working directory acts as a hidden parameter.
        """
        if wk_dir is None:
            self._wk_dir = os.getcwd()
        else:
            self._wk_dir = wk_dir

    @property
    def temp_file(self):
        return self._temp_file

    @temp_file.setter
    def temp_file(self, temp_file):
        self._temp_file = temp_file

    @property
    def tree(self):
        return self._tree

    @tree.setter
    def tree(self, tree):
        self._tree = tree

    @property
    def sim_map_name(self):
        return self._sim_map_name

    @sim_map_name.setter
    def sim_map_name(self, sim_map_name):
        try:
            self._sim_map_name.update(sim_map_name)
        except AttributeError:
            self._sim_map_name = sim_map_name

    @property
    def sim_map(self):
        return self._sim_map

    @sim_map.setter
    def sim_map(self, sim_map):
        try:
            self._sim_map.update(sim_map)
        except AttributeError:
            self._sim_map = sim_map
    
    def __call__(self, raiseonerr=True):
        """
        The main routine for running hammurabiX executable.
        
        Parameters
        ----------
        
        raiseonerr : bool
            Error raising flag, if not raising on error, do logging.
        """
        # create new temp parameter file
        if self.temp_file is self._base_file:
            self._new_xml_copy()
        if raiseonerr:
            logfile = subprocess.PIPE
            errfile = subprocess.PIPE
        else:
            logfile = open('run.log', 'w')
            errfile = open('err.log', 'w')
        try:
            temp_process = subprocess.Popen([self._executable, self._temp_file],
                                            stdout=logfile,
                                            stderr=errfile)
            returncode = temp_process.wait()
        except Exception as e:
            if raiseonerr:
                raise OSError("ERROR subprocess.Popen raised {}".format(e))
            else:
                print("WARNING subprocess.Popen raised {}".format(e))
        if returncode != 0:
            # This only returns the output if you used subprocess.PIPE, i.e. writelog=False
            last_call_log, last_call_err = temp_process.communicate()
            print(last_call_log)
            print(last_call_err)
            if raiseonerr:
                raise ValueError("HamX returned status {}".format(returncode))
        if not raiseonerr:
            logfile.close()
            errfile.close()
        # grab output maps and delete temp files
        self._get_sims()
        self._del_xml_copy()

    def _new_xml_copy(self):
        """
        Make a temporary parameter file copy and rename output file with random mark.
        """
        # create a random file name which doesn't exist currently
        fd, new_path = tempfile.mkstemp(prefix='params_', suffix='.xml', dir=self._wk_dir)
        os.close(fd)
        # retrieve the temporary mark
        rnd_idx = new_path[new_path.index('params_')+7:-4]
        self.temp_file = new_path
        # copy base_file to temp_file
        root = self.tree.getroot()
        # count number of sync output files
        if root.find("./observable/sync[@cue='1']") is not None:
            self._do_sync = True
            # for each sync output
            for sync in root.findall("./observable/sync[@cue='1']"):
                freq = str(sync.get('freq'))
                nside = str(sync.get('nside'))
                self.sim_map_name[('sync', freq, nside)] = os.path.join(self.wk_dir,
                                                                        'sync_'+freq+'_'+nside+'_'+rnd_idx+'.bin')
                sync.set('filename', self.sim_map_name[('sync', freq, nside)])
        fd = root.find("./observable/faraday[@cue='1']")
        if fd is not None:
            self._do_fd = True
            nside = str(fd.get('nside'))
            self.sim_map_name[('fd', 'nan', nside)] = os.path.join(self.wk_dir, 'fd_'+nside+'_'+rnd_idx+'.bin')
            fd.set('filename', self.sim_map_name[('fd', 'nan', nside)])
        dm = root.find("./observable/dm[@cue='1']")
        if dm is not None:
            self._do_dm = True
            nside = str(dm.get('nside'))
            self.sim_map_name[('dm', 'nan', nside)] = os.path.join(self._wk_dir, 'dm_'+nside+'_'+rnd_idx+'.bin')
            dm.set('filename', self.sim_map_name[('dm', 'nan', nside)])
        # automatically create a new file
        self.tree.write(self.temp_file)

    def _get_sims(self):
        """
        Grab simulation output from disk and delete corresponding fits file.
        """
        # locate the keys
        sync_key = list()
        dm_key = None
        fd_key = None
        for k in self.sim_map_name.keys():
            if k[0] == 'dm':
                dm_key = k
            elif k[0] == 'fd':
                fd_key = k
            elif k[0] == 'sync':
                sync_key.append(k)
            else:
                raise ValueError('mismatched key %s' % str(k))
        # read dispersion measure and delete file
        if self._do_dm is True:
            self.sim_map[(dm_key[0], dm_key[1], dm_key[2], 'nan')] = self._read_del(self.sim_map_name[dm_key])
        # read faraday depth and delete file
        if self._do_fd is True:
            self.sim_map[(fd_key[0], fd_key[1], fd_key[2], 'nan')] = self._read_del(self.sim_map_name[fd_key])
        # read synchrotron pol. and delete file
        if self._do_sync is True:
            for i in sync_key:
                self.sim_map[(i[0], i[1], i[2], 'I')] = self._read_del(self.sim_map_name[i][:-4]+'_I.bin')
                self.sim_map[(i[0], i[1], i[2], 'Q')] = self._read_del(self.sim_map_name[i][:-4]+'_Q.bin')
                self.sim_map[(i[0], i[1], i[2], 'U')] = self._read_del(self.sim_map_name[i][:-4]+'_U.bin')
                # polarisation intensity
                #self.sim_map[(i[0], i[1], i[2], 'PI')] = np.sqrt(np.square(Qs) + np.square(Us))
                # polarisatioin angle, IAU convention
                #self.sim_map[(i[0], i[1], i[2], 'PA')] = np.arctan2(Us, Qs)/2.0
                    
    def _read_del(self, path):
        """
        Read a single binary file into a numpy array,
        and then delete the file.
        """
        assert(os.path.isfile(path))
        loaded_map = np.fromfile(path, dtype=np.float64)
        os.remove(path)
        return loaded_map

    def _del_xml_copy(self):
        """
        Delete temporary parameter file copy.
        """
        if self.temp_file is self._base_file:
            raise ValueError('read only')
        else:
            os.remove(self._temp_file)
            self._temp_file = self._base_file

    """
    THE FOLLOWING FUNCTIONS ARE RELATED TO XML FILE MANIPULATION
    """

    def mod_par(self, keychain=None, attrib=None):
        """
        Modify parameter in the XML tree,
        argument of type ['path','to','target'], {attrib}.
        
        Parameters
        ----------
        
        keychain : list
            A list of key entry strings.
        
        attrib : dict
            A dict in format {'tag': 'content'},
            where 'tag' is know as XML attribute name, of given XML element,
            followed by the 'content' which records the attribute value.
            If attribute 'tag' already exists, then new attrib will be assigned.
            If attribute 'tag' is not found, then new 'tag' will be inserted.
        """
        # input type check
        if type(attrib) is not dict or type(keychain) is not list:
            raise ValueError('wrong input %s %s' % (keychain, attrib))
        root = self.tree.getroot()
        path_str = '.'
        for key in keychain:
            path_str += '/' + key
        target = root.find(path_str)
        if target is None:
            raise ValueError('wrong path %s' % path_str)
        for i in attrib:
            target.set(i, attrib.get(i))

    def add_par(self, keychain=None, subkey=None, attrib=None):
        """
        Add new subkey under keychain in the XML tree,
        argument of type ['path','to','target'], 'subkey', {attrib},
        or of type ['path','to','target'], 'subkey'.
        """
        # input type check
        if type(keychain) is not list or type(subkey) is not str:
            raise ValueError('wrong input %s %s %s' % (keychain, subkey, attrib))
        if attrib is not None and type(attrib) is dict:
            root = self.tree.getroot()
            path_str = '.'
            for key in keychain:
                path_str += '/' + key
            target = root.find(path_str)
            # check if (subkey, attrib) already exists
            for existed in target.findall(subkey):
                if existed.attrib == attrib:
                    raise ValueError('repeatitive input %s %s %s' % (keychain, subkey, attrib))
            # if not, add to the given position
            et.SubElement(target, subkey, attrib)
        elif attrib is None:
            root = self.tree.getroot()
            path_str = '.'
            for key in keychain:
                path_str += '/' + key
            target = root.find(path_str)
            et.SubElement(target, subkey)
        else:
            raise ValueError('wrong input %s %s %s' % (keychain, subkey, attrib))

    def print_par(self, keychain=None):
        """
        Print a certain parameter,
        argument of type ['path','to','key'] (e.g. ['grid','observer','x']),
        print all parameters down to the keychain children level.
        """
        # input type check
        if type(keychain) is not list:
            raise ValueError('wrong input %s' % keychain)
        root = self.tree.getroot()
        # print top parameter level if no input is given
        if keychain is None:
            for child in root:
                print(child.tag, child.attrib)
        else:
            path_str = '.' 
            for key in keychain:
                path_str += '/' + key
            for target in root.findall(path_str):
                print(target.tag, target.attrib)
                for child in target:
                    print('|--> ', child.tag, child.attrib)

    def del_par(self, keychain=None, clean=False):
        """
        Deletes an parameter and all of its children,
        argument of type ['keys','to','target'] (e.g. ['grid','observer','x']),
        if clean=True, delete all parameters that match given keychain.
        """
        # input type check
        if type(keychain) is not list:
            raise ValueError('wrong input %s' % keychain)
        assert isinstance(clean, bool)
        root = self.tree.getroot()
        if keychain is not None:
            path_str = '.'
            par_path_str = '.'
            n = 1
            for key in keychain:
                path_str += '/' + key
                n += 1
                if n is len(keychain):
                    par_path_str = path_str
            target = root.find(path_str)
            parent = root.find(par_path_str)
            if target is None or parent is None:
                raise ValueError('wrong path %s' % path_str)
            if clean:
                for i in root.findall(path_str):
                    parent.remove(i)
            else:
                parent.remove(target)
        else:
            raise ValueError('empty keychain')
