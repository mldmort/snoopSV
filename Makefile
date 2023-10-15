
install:
	python setup.py install

test:
	cd tests/unit
	pytest
