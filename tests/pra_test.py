import unittest

import modular_core.libfundamental as lfu
lfu.USING_GUI = False
#import modular_pra as pra

import os, sys, pdb

#log = open(os.path.join(os.getcwd(), 'test_ensemble.log'), 'w')
#sys.stdout = log

class dummyTestCase(unittest.TestCase):
	#simple_mcfg = os.path.join(os.getcwd(), 
	#			'stringchemical_dep_mcfgs', 
	#			'MM_kinetics_boring.mcfg')

	def pause(self, *args, **kwargs):
		sys.stdout = sys.__stdout__
		pdb.set_trace()
		sys.stdout = log

	def test_can_import(self):
                import modular_pra as pra
		self.assertFalse(pra == None)

if __name__ == '__main__':
    unittest.main()


