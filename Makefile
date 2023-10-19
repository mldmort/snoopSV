
install:
	python -m pip install -r requirements.txt .

dev:
	python -m pip install -e .

test:
	cd tests/unit
	pytest
