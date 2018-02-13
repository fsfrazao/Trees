from setuptools import setup

setup(name='trees_ibm',
      version='0.2.3',
      description='Individual-based tree community models',
      url='http://github.com/fsfrazao/Trees',
      author='Fabio Frazao',
      author_email='fsfrazao@gmail.com',
      license='GNU General Public License v3.0',
      packages=['trees_ibm'],
      install_requires=[
          'numpy',
          'py_ibm',
          ],
      setup_requires=['pytest-runner',],
      tests_require=['pytest',],
      zip_safe=False)
