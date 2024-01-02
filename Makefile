# Make variables
#----------------------------------------------------------------------
BOLD_CYAN = \033[1;36m
BOLD_GREEN := \033[1;32m
GREEN = \033[0;32m
NO_COLOR = \033[0m


## General
# ----------------------------------------------------------------------
help:  ## Print this help message
	@egrep -h '(\s|^)##\s' $(MAKEFILE_LIST) \
	| sed -E "s/^## (.*)/\n$$(printf "${BOLD_CYAN}")\1$$(printf "${NO_COLOR}")/g" \
	| awk 'BEGIN {FS = ":.*?## "}; {printf "${GREEN}%-25s${NO_COLOR} %s\n", $$1, $$2}'
	@echo


## Development environment
# ----------------------------------------------------------------------
PACKAGE_NAME = yhaplo
ENV_NAME = $(PACKAGE_NAME)

dev-pyenv-virtualenv:  ## Set up pyenv-virtual-env-based development environment
	pyenv uninstall --force $(ENV_NAME)
	pyenv local --unset
	pyenv virtualenv $(ENV_NAME)
	pyenv local $(ENV_NAME)
	pip install --upgrade pip setuptools wheel
	$(MAKE) dev-install
	$(MAKE) init-hooks

dev-install:  ## Install package as editable, with all optional dependencies
	pip install --editable .[dev]
	python -m ipykernel install --user --name $(ENV_NAME) --display-name $(PACKAGE_NAME)


## Pre-commit hooks
# ----------------------------------------------------------------------
init-hooks: install-hooks update-hooks  ## Install and update hooks

install-hooks:  ## Install hooks
	pre-commit install --install-hooks --hook-type pre-commit --hook-type commit-msg

update-hooks:  ## Update hooks
	pre-commit autoupdate

run-hooks:  ## Run hooks
	pre-commit run

run-hooks-all:  ## Run hooks on all files
	pre-commit run --all-files

lint: run-hooks-all  ## Alias for run-hooks-all


## Testing
# ----------------------------------------------------------------------
test:  ## Run unit tests
	pytest


