# Makefile for running tests on the python code here

# These variables can be overridden from the command-line
python = not-set
verbose = not-set
debug = not-set

ifneq ($(python), not-set)
    PYTHON=$(python)
else
    PYTHON=python3
endif

ifneq ($(debug), not-set)
    TEST_ARGS+=--debug
endif
ifneq ($(verbose), not-set)
    TEST_ARGS+=--verbose
endif

PYLINT=pylint
PYLINT_ARGS=-j 4 --rcfile=fates/.pylintrc
PYLINT_SRC = \
	fates ../testing

all: test black lint
# ----------------------------------------------------------------
# See the stest target about this issue
	@echo "Run './run_fates_py_tests --sys' by hand afterwards"
# ----------------------------------------------------------------
	@echo
	@echo
	@echo "Successfully ran all standard tests"

test: utest stest

.PHONY: utest
utest: FORCE
	$(PYTHON) ./run_fates_py_tests $(TEST_ARGS) --unit

.PHONY: stest
stest: FORCE
# ----------------------------------------------------------------
# EBK 2024-03-19: Comment out running here because of this issue:
# https://github.com/ESCOMP/fates/pull/2363#issuecomment-1967884908
#$(PYTHON) ./run_fates_py_tests $(TEST_ARGS) --sys
# Instead run by hand which seems to be working for now...
# ----------------------------------------------------------------
	@echo "System tests currently don't run under Make so..."
	@echo "Run './run_fates_py_tests --sys' by hand afterwards"

.PHONY: lint
lint: FORCE
	$(PYLINT) $(PYLINT_ARGS) $(PYLINT_SRC)

.PHONY: black
# Run the black check on all of the python files here and undeneath.
# Use the black configure file to explicitly set a few things and specify the exact files.
black: FORCE
	black --check --config pyproject.toml . ../testing

.PHONY: run_black
# Run black on all of the python files here and undeneath.
# Use the black configure file to explicitly set a few things and specify the exact files.
run_black: FORCE
	black --config pyproject.toml . ../testing

.PHONY: clean
clean: FORCE
	find . -name '*.pyc' -exec rm {} \;

FORCE:

