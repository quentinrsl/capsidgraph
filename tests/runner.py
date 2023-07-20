import unittest
import test_analyser
import test_generation

suite = unittest.TestSuite()
loader = unittest.TestLoader()

suite.addTests(loader.loadTestsFromModule(test_analyser))
suite.addTests(loader.loadTestsFromModule(test_generation))

runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)