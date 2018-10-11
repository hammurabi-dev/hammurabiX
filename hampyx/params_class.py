# -*- coding: utf-8 -*-
'''
Contributed by Joe Taylor
'''

import os
import subprocess
import healpy as hp
import xml.etree.ElementTree as et
import random as rnd
import numpy as np
import tempfile as tf
from shutil import copyfile

class Params(object):
		
	# we try to keep __init__ as light as possible
	def __init__(self,input_file=None):
		
		# for simplicity we nail down the working directory
		#self.dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
		self.dir = os.path.abspath('.')
		self.read(input_file)


	def read(self,parameter_file=None):
		import sys
		import time	
		if type(parameter_file) is str:
			# nail down parameter file
			self.param_file = os.path.join(self.dir,parameter_file)
		else:
			if parameter_file is not None:
				print ('Error: Invalid parameter file specified, returning default parameters')
			self.param_file = os.path.join(self.dir,'defaults.xml')

		self.tree = et.parse(self.param_file)
		

	def write(self,filename=None):
	
		if type(filename) is str:
			self.tree.write(filename)
		else:
			print ('Error: Please input a string for the filename to write to')


	def add_ele(self, keychain): 	#keychain has the form [['path','to','target'],'tag','attrib']
		
		# creates an element as a subelement of a parent with a tag and attributes

		root = self.tree.getroot()
		parentkeys = keychain[0]
		eletag = keychain[1]
		attrib = keychain[2]
						
		if parentkeys is not None and type(eletag) is str:
			path_str = '.'
			for key in parentkeys:
				path_str += '/' + key
			target = root.find(path_str)
			if type(attrib) is dict:
				et.SubElement(target,eletag,attrib)
			else:
				et.SubElement(target,eletag)

		else:
			print ('Error: Valid input required')


	def del_ele(self, keys=None):

		# deletes an element and all of its children
		# keys is the path to the element (e.g. ['Grid','SunPosition','x'])

		root = self.tree.getroot()

		if keys is not None:
			path_str = '.'
			n = 1
			for key in keys:
				if n is len(keys):
					par_path_str = path_str
				path_str += '/' + key
				n += 1
			target = root.find(path_str)
			parent = root.find(par_path_str)
			parent.remove(target)

		else:
			print ('Error: Input required')

	def get_ele(self, keys=None, opts=None):
                
                root = self.tree.getroot()
            
                # print top parameter level if no input is given
                if keys is None:
                	for child in root:
                        	print child.tag, child.attrib
                
                else:

                    if opts is None: # print the selected level
                    
                        path_str = '.' 
    		        for key in keys:
		        	path_str += '/' + key
    		
		        target = root.find(path_str)
                        print target.tag, target.attrib
                        for child in target:
                            print '|-->', child.tag, child.attrib

                    elif opts is 'All': # print all levels down to the selected level

                        path_str = '.' 
    		        for key in keys:
			    path_str += '/' + key
                            target = root.find(path_str)
                            print target.tag
                            for child in target:
                                print ('|-->', child.tag, child.attrib)                        
                        
                    else:
                        print ('Error: option is not supported')

 
        #def mod_par(self, keys=None, tag=None, attrib=None):
        def mod_par(self,keylist=None):

		if len(keylist) is 2:
			keys = keylist[0]
			tag = None
			attrib = keylist[1]
		elif len(keylist) is 3:
			keys = keylist[0]
			tag = keylist[1]
			attrib = keylist[2]
		else:
			print ('Error: invalid input')
			exit(1)

                if (keys is None or attrib is None):
                    print ('Error: missing input')
                    exit(1)
                
                if tag is None:
                    tag = 'value'
		
		root = self.tree.getroot()
		
		path_str = '.'
		for key in keys:
			path_str += '/' + key
    		
		target = root.find(path_str)
		if target is None:
			print ('Error: wrong element path:')
			print path_str
			exit(1)
		
		target.set(tag,attrib)


	def mod_all(self,phrase=None,attribname=None,newval=None):
		
		if phrase is None:
			print ('Error: no input given')	
		elif type(phrase) is str:
			root = self.tree.getroot()
			for target in root.iter(phrase):
				if newval is not None:
					target.set(attribname,newval)
		else:
			print ('Error: invalid input given')
