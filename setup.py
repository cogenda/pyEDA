from distutils.core import setup

setup(name='Distutils',
      version='1.0',
      description='EDA/TCAD framework and utilities for education purposes.',
      author='Shen Chen',
      author_email='shenchen@cogenda.com',
      url='http://www.cogenda.com',
      packages=['pyEDA', 'pyEDA.Circuit', 'pyEDA.Compact', 'pyEDA.FVMEqn', 'pyEDA.Mesh', 'pyEDA.PDE', 'pyEDA.Device']
     )
