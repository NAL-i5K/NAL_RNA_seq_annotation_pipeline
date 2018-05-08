import unittest
import subprocess


class JavaSetupTestCase(unittest.TestCase):
    def test(self):
        subprocess.run(['java', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)


class SRAToolkitSetupTestCase(unittest.TestCase):
    def test(self):
        subprocess.run(['fastq-dump', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

