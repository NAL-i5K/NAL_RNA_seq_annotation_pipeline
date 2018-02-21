import unittest
import subprocess


class SetupTestCase(unittest.TestCase):
    def test(self):
        subprocess.run(['java', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)