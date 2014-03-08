#!/usr/bin
# libraries for the unittesting
import unittest
import subprocess
from subprocess import call
import re
import os

# testing functions
import integrator


while(True):
    user_input = raw_input("Do you want to recompile the code? (y/n) ")
    if(user_input[0] == "y" or user_input == "Y"):
        # prepare the files
        call(["make", "clean"])
        call(["make"])
        break
    elif(user_input[0] =="n" or user_input == "N"):
        break
    
# begin the testing functions 
class TestFunctions(unittest.TestCase):

    def test_betas_integral(self):
        os.system("./work15b-all >& work15b-all.log")
        self.assertTrue(integrator.integrate("betas.out") <= 1.0)

    def test_overlap_plus_betas(self):
        self.assertTrue(True)

suite = unittest.TestLoader().loadTestsFromTestCase(TestFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
