import unittest
from atlasdavinci.cli import main

class TestCLI(unittest.TestCase):
    def test_main(self):
        """
        Test the main function
        """
        # 这里我们只是确保main函数可以运行而不抛出异常
        try:
            main()
        except Exception as e:
            self.fail(f"main() raised {type(e).__name__} unexpectedly!")

if __name__ == '__main__':
    unittest.main() 