import subprocess
import sys

def install(package):
    subprocess.call([sys.executable, "-m", "pip", "install", package])

dependencies = ['numpy', 'pandas']

# Example
if __name__ == '__main__':
	for d in dependencies:
		install(d)