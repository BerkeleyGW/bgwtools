from distutils.core import setup

setup(
	name='BGWTools',
	version='0.9.0',
	author='Felipe Homrich da Jornada',
	author_email='jornada@civet.berkeley.edu',
	packages=['towelstuff', 'towelstuff.test'],
#	scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
#	url='http://pypi.python.org/pypi/TowelStuff/',
	license='LICENSE.txt',
	description='BerkeleyGW Tools and Scripts',
	long_description=open('README.txt').read(),
	install_requires=[
	],
)
