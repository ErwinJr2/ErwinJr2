import os
import sys

os.chdir(os.path.dirname(__file__))

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

if __name__ == "__main__":
    print(os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..')))
