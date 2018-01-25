from setuptools import setup

setup(name='plot_nanopore_run',
      version='0.1',
      description='Generates plots reporting on nanopore runs',
      url='https://github.com/jopineda/plot-nanopore-run',
      author='Joanna Pineda',
      author_email='joanna.pineda@mail.mcgill.ca',
      license='MIT',
      python_requires='>=2.7',
      packages=['plot_nanopore_run'],
      setup_requires=['numpy'],
      install_requires=[
            'gevent',
            'greenlet',
            'matplotlib==2.0.0'
      ],
      zip_safe=False)
