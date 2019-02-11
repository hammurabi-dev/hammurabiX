# -*- coding: utf-8 -*-
'''
Developed by Joe Taylor, based on the initial work by Theo Steininger
Reviewed by Jiaxin Wang

temporary/testing version of hammurabiX python wrapper

warning:
working directory is set as the same directory as this file
it relies on subprocess to fork c++ routine and passing data through disk
so, it is not fast

methods:

# Import class
In []: import hampyx as hpx

# Initialize object
In []: object = hpx.hampyx (exe_path, xml_path)

hammurabiX executable path is by default '/usr/local/hammurabi/bin/hamx'
while xml file path is by default './'

# Modify parameter value from base xml file to temp xml file
In []: object.mod_par (keychain=['key1','key2',...], attrib={'tag':'content'})

# Add new parameter with or without attributes
In []: object.add_par (keychain=['key1','key2',...], subkey='keyfinal', attrib={'tag':'content'})

the new parameter subkey, will be added at the path defined by keychain

# Delete parameter
In []: object.del_par (keychain=['key1','key2',...])

if additional argument opt='all', then all matching parameters will be deleted

the strings 'key1', 'key2', etc represent the path to the desired parameter, going through the xml
the "tag" is the label for the parameter: eg. "Value" or "cue" or "type"
the "content" is the content under the tag: eg. the string for the tag "filename"

# Look through the parameter tree in python
In []: object.print_par(keychain=['key1','key2',...])

this will return the current value of the parameter in the XML associated with the path "key1/key2/.../keyfinal/"

# Run the executable
In []: object.call()

if additional verbose=True (by default is False) hampyx_run.log and hampyx_err.log will be dumped to disk
notice that dumping logs is not thread safe, use quiet mode in threading

after this main routine, object.sim_map will be filled with simulation outputs from hammurabiX
the structure of object.sim_map contains arrays under entries:

object.sim_map['sync']['frequency']['I'] # synchrotron intensity map at 'frequency' 
object.sim_map['sync']['frequency']['Q'] # synchrotron Q map at 'frequency' 
object.sim_map['sync']['frequency']['U'] # synchrotron U map at 'frequency' 
object.sim_map['sync']['frequency']['PI'] # synchrotron pol. intensity at 'frequency' 
object.sim_map['sync']['frequency']['PA'] # synchrotron pol. angle at 'frequency' (IAU convention)
object.sim_map['fd'] # Faraday depth map
object.sim_map['dm'] # dispersion measure map

detailed caption of each function can be found with their implementation
'''

import os
import sys
import time
import subprocess
import healpy as hp
import xml.etree.ElementTree as et
import random as rnd
import numpy as np
import tempfile as tf

class hampyx(object):
    
    '''
    default executable path is '/usr/local/hammurabi/bin/hamx'
    default executable path is './params.xml'
    '''
    def __init__(self, exe_path=None, xml_path=None):
        # current working directory
        self.wk_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        # encapsulated below
        self.exe_path = exe_path
        self.xml_path = xml_path
        # assign tmp file name by base file name
        self.temp_file = self.base_file
        # read from base parameter file
        self.tree = et.parse(self.base_file)
        # simulation output
        self.sim_map_name = {}
        self.sim_map = {}
        # 'sync' is nested dict
        self.sim_map_name['sync'] = {}
        self.sim_map['sync'] = {}
        # switches
        self.do_sync = False
        self.do_dm = False
        self.do_fd = False
    
    @property
    def exe_path(self):
        return self._exe_path
    
    @property
    def xml_path(self):
        return self._xml_path
    
    '''
    by default hammurabiX executable "hamx" is in the same directory as this file
    '''
    @exe_path.setter
    def exe_path(self,exe_path):
        if exe_path is None:#{
            self._exe_path = os.path.abspath('/usr/local/hammurabi/bin/hamx')
            self.executable = self.exe_path
        #}
        elif isinstance(exe_path, basestring):#{
            self._exe_path = exe_path
            self.executable = os.path.abspath(self.exe_path)
        #}
        else:
            raise TypeError (exe_path + ' is expected to be of basestring type')

    @xml_path.setter
    def xml_path(self,xml_path):
        if xml_path is None:#{
            self._xml_path = os.path.join(self.wk_dir,'params.xml')
            self.base_file = self.xml_path
        #}
        elif isinstance(xml_path, basestring):#{
            self._xml_path = os.path.abspath(xml_path)
            self.base_file = self.xml_path
        #}
        else:
            raise TypeError (xml_path + ' is expected to be of basestring type')

    '''
    the main routine for running hammurabiX executable
    '''
    def __call__(self, verbose=False):
        # create new temp parameter file
        if self.temp_file is self.base_file:#{
            self._new_xml_copy()
        #}
        # if need verbose output
        if verbose is True:#{
            logfile = open('hampyx_run.log', 'w')
            errfile = open('hampyx_err.log','w')
            temp_process = subprocess.Popen([self.executable,self.temp_file],stdout=logfile,stderr=errfile)
            temp_process.wait()
            logfile.close()
            errfile.close()
        #}
        # if quiet, only print upon error
        else:#{
            temp_process = subprocess.Popen([self.executable,self.temp_file],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            temp_process.wait()
            if temp_process.returncode != 0:#{
                last_call_log,last_call_err = temp_process.communicate()
                print last_call_log
                print last_call_err
            #}
        #}
        # grab output maps and delete temp files
        self._get_sims()
        self._del_xml_copy()

    '''
    grab simulation output from disk and delete corresponding fits file
    '''
    def _get_sims(self):
        # locate the files
        dm_path = self.sim_map_name['dm']
        fd_path = self.sim_map_name['fd']
        sync_paths = {}
        # sync_path has multiple contents
        for i in self.sim_map_name['sync']:#{
            sync_paths[i] = self.sim_map_name['sync'].get(i)
        #}
        # read dispersion measure and delete file
        if self.do_dm is True:#{
            if(os.path.isfile(dm_path)):#{
                [DM] = self._read_fits_file(dm_path)
                self.sim_map['dm'] = DM
                os.remove(dm_path)
            #}
            else:#{
                raise ValueError (dm_path + ' not found')
            #}
        #}
        # read faraday depth and delete file
        if self.do_fd is True:#{
            if(os.path.isfile(fd_path)):#{
                [Fd] = self._read_fits_file(fd_path)
                self.sim_map['fd'] = Fd
                os.remove(fd_path)
            #}
            else:#{
                raise ValueError (fd_path + ' not found')
            #}
        #}
        # read synchrotron pol. and delete file
        if self.do_sync is True:#{
            for i in sync_paths:#{
                # if file exists
                if(os.path.isfile(sync_paths.get(i))):#{
                    [Is,Qs,Us] = self._read_fits_file(sync_paths.get(i))
                    # build top level map for nesting
                    self.sim_map['sync'][i] = {}
                    self.sim_map['sync'][i]['I'] = Is
                    self.sim_map['sync'][i]['Q'] = Qs
                    self.sim_map['sync'][i]['U'] = Us
                    # polarisation intensity
                    self.sim_map['sync'][i]['PI'] = np.sqrt(np.square(Qs) + np.square(Us))
                    # polarisatioin angle, IAU convention
                    self.sim_map['sync'][i]['PA'] = np.arctan2(Us,Qs)/2.0
                    os.remove(sync_paths.get(i))
                #}
                else:#{
                    raise ValueError (sync_paths.get(i) + ' not found')
                #}
            #}
        #}

    '''
    read a single fits file with healpy
    '''
    def _read_fits_file(self, path):
        rslt = []
        i = 0
        while True:#{
            try:#{
                loaded_map = hp.read_map(path,verbose=False,field=i)
                rslt += [loaded_map]
                i += 1
            #}
            except IndexError:#{
                break
            #}
        #}
        return rslt

    '''
    make a temporary parameter file copy and rename output file with random mark
    '''
    def _new_xml_copy(self):
        # create a random file name which doesn't exist currently
        fd,new_path = tf.mkstemp(prefix='params_',suffix='.xml',dir=self.wk_dir)
        os.close(fd)
        rnd_idx = new_path[new_path.index('params_')+7:-4]
        self.temp_file = new_path
        # copy base_file to temp_file
        root = self.tree.getroot()
        # count number of sync output files
        if root.find("./Obsout/Sync[@cue='1']") is not None:#{
            self.do_sync = True
            # for each sync output
            for sync in root.findall("./Obsout/Sync[@cue='1']"):#{
                freq = str(sync.get('freq'))
                self.sim_map_name['sync'][freq] = os.path.join(self.wk_dir,'iqu_sync_'+freq+'_'+rnd_idx+'.fits')
                sync.set('filename',self.sim_map_name['sync'][freq])
            #}
        #}
        self.sim_map_name['fd'] = os.path.join(self.wk_dir,'fd_'+rnd_idx+'.fits')
        # in case no fd map is required
        if root.find("./Obsout/Faraday[@cue='1']") is not None:#{
            self.do_fd = True
            root.find("./Obsout/Faraday[@cue='1']").set('filename',self.sim_map_name['fd'])
        #}
        self.sim_map_name['dm'] = os.path.join(self.wk_dir,'dm_'+rnd_idx+'.fits')
        # in case no dm map is required
        if root.find("./Obsout/DM[@cue='1']") is not None:#{
            self.do_dm = True
            root.find('./Obsout/DM').set('filename',self.sim_map_name['dm'])
        #}
        # automatically create a new file
        self.tree.write(self.temp_file)
	
    '''
    delete temporary parameter file copy
    '''
    def _del_xml_copy(self):
        if self.temp_file is self.base_file:#{
            raise ValueError (self.temp_file + ' read only')
        #}
        else:#{
            os.remove(self.temp_file)
            self.temp_file = self.base_file
        #}

    '''
    THE FOLLOWING FUNCTIONS ARE RELATED TO XML FILE MANIPULATION
    '''

    '''
    modify parameter in self.tree
    argument of type ['path','to','target'], {attrib}
    attrib of type {'tag': 'content'}
    if attribute 'tag' already exists, then new attrib will be assigned
    if attribute 'tag' is not found, then new 'tag' will be inserted
    '''
    def mod_par(self, keychain=None, attrib=None):
        # input type check
        if type(attrib) is not dict or type(keychain) is not list:#{
            raise ValueError ('wrong input ',keychain,attrib)
        #}
        root = self.tree.getroot()
        path_str = '.'
        for key in keychain:#{
            path_str += '/' + key
        #}
        target = root.find(path_str)
        if target is None:#{
            raise ValueError ('wrong path ',path_str)
        #}
        for i in attrib:#{
            target.set(i,attrib.get(i))
        #}

    '''
    add new subkey under keychain in the tree
    argument of type ['path','to','target'], 'subkey', {attrib}
    or of type ['path','to','target'], 'subkey'
    '''
    def add_par(self, keychain=None, subkey=None, attrib=None):
        # input type check
        if type(keychain) is not list or type(subkey) is not str:#{
            raise ValueError ('wrong input ',keychain,subkey,attrib)
        #}
        if attrib is not None and type(attrib) is dict:#{
            root = self.tree.getroot()
            path_str = '.'
            for key in keychain:#{
                path_str += '/' + key
            #}
            target = root.find(path_str)
            et.SubElement(target,subkey,attrib)
        #}
        elif attrib is None:#{
            root = self.tree.getroot()
            path_str = '.'
            for key in keychain:#{
                path_str += '/' + key
            #}
            target = root.find(path_str)
            et.SubElement(target,subkey)
        #}
        else:#{
            raise ValueError ('wrong input ',keychain,subkey,attrib)
        #}

    '''        
    print a certain parameter
    argument of type ['path','to','key'] (e.g. ['Grid','SunPosition','x'])
    print all parameters down to the keychain children level
    '''
    def print_par(self, keychain=None):
        # input type check
        if type(keychain) is not list:#{
            raise ValueError ('wrong input ',keychain)
        #}
        root = self.tree.getroot()
        # print top parameter level if no input is given
        if keychain is None:#{
            for child in root:#{
                print child.tag, child.attrib
            #}
        #}
        else:#{
            path_str = '.' 
            for key in keychain:#{
                path_str += '/' + key
            #}
            #target = root.find(path_str)
            for target in root.findall(path_str):#{
                print target.tag, target.attrib
                for child in target:#{
                    print ('|--> ', child.tag, child.attrib)
            #}
        #}
    #}

    '''        
    deletes an parameter and all of its children
    argument of type ['keys','to','target'] (e.g. ['Grid','SunPosition','x'])
    if opt='all', delete all parameters that match given keychain
    '''
    def del_par(self, keychain=None, opt=None):
        # input type check
        if type(keychain) is not list:#{
            raise ValueError ('wrong input ',keychain,opt)
        #}
        root = self.tree.getroot()
        if keychain is not None:#{
            path_str = '.'
            n = 1
            for key in keychain:#{
                path_str += '/' + key
                n += 1
                if n is len(keychain):#{
                    par_path_str = path_str
                #}
            #}
            target = root.find(path_str)
            parent = root.find(par_path_str)
            if target is None or parent is None:#{
                raise ValueError ('wrong path ',path_str)
            #}
            if opt is None:#{
                parent.remove(target)
            #}
            elif opt is 'all':#{
                for i in root.findall(path_str):#{
                    parent.remove(i)
                #}
            #}
            else:#{
                raise ValueError ('unsupported option ',keychain,opt)
            #}
        #}
        else:#{
            raise ValueError ('empty keychain ')
        #}
