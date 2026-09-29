PYTHON ?= python3
export PYTHONPATH := .

.PHONY: all analyze test check clean

all: analyze test

analyze:
	$(PYTHON) scripts/run_verification.py

test:
	$(PYTHON) -m unittest discover -s tests -v

check:
	$(PYTHON) -c "import sys, numpy; print(sys.version); print('numpy', numpy.__version__)"

clean:
	$(RM) results/*.csv results/*.json figures/*.svg

