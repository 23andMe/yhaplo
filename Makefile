# Make variables
#----------------------------------------------------------------------
CYAN := \033[0;36m
GREEN := \033[0;32m
BOLD_CYAN := \033[1;36m
BOLD_GREEN := \033[1;32m
RESET_COLOR := \033[0m


## General
# ----------------------------------------------------------------------
help:  ## Print this help message
	@egrep -h '(\s|^)##\s' $(MAKEFILE_LIST) \
	| sed -E "s/^## (.*)/\n$$(printf "${BOLD_GREEN}")\1$$(printf "${RESET_COLOR}")/g" \
	| awk 'BEGIN {FS = ":.*?## "}; {printf "${CYAN}%-25s${RESET_COLOR} %s\n", $$1, $$2}'
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
	pip install --upgrade uv
	uv pip install --python=$(shell which python) --upgrade pip setuptools wheel
	$(MAKE) dev-install
	$(MAKE) dev-jupyter
	$(MAKE) init-hooks

dev-install:  ## Install package as editable, with all optional dependencies
	@printf "\n${BOLD_GREEN}Installing package as editable, with all optional dependencies${RESET_COLOR}...\n\n"
	uv pip install --python=$(shell which python) --editable .[dev]

dev-jupyter:  ## Add Jupyter kernel
	@printf "\n${BOLD_GREEN}Installing Jupyter kernel${RESET_COLOR}...\n\n"
	python -m ipykernel install --user --name $(ENV_NAME) --display-name $(PACKAGE_NAME)


## Pre-commit hooks
# ----------------------------------------------------------------------
init-hooks: install-hooks update-hooks  ## Install and update hooks

install-hooks:  ## Install hooks
	@printf "\n${BOLD_GREEN}Installing pre-commit hooks${RESET_COLOR}...\n\n"
	pre-commit install --install-hooks --hook-type pre-commit --hook-type commit-msg

update-hooks:  ## Update hooks
	@printf "\n${BOLD_GREEN}Upgrading pre-commit hooks${RESET_COLOR}...\n\n"
	pre-commit autoupdate

run-hooks:  ## Run hooks
	pre-commit run

run-hooks-all:  ## Run hooks on all files
	pre-commit run --all-files


## Linting
# ----------------------------------------------------------------------
lint: run-hooks-all  ## Alias for run-hooks-all

ruff:  ## Run Ruff linting and formatting, fixing violations
	ruff check --fix --unsafe-fixes
	ruff format


## Testing
# ----------------------------------------------------------------------
test:  ## Run unit tests
	pytest


