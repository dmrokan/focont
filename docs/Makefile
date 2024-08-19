doc_path := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
project_path := $(doc_path)/..
sphinx_path := $(project_path)/_docs

doc:
	cd "$(doc_path)/.." && python3 -m pip install -e '.[doc]'
	sphinx-apidoc -o "$(sphinx_path)/source/" "$(project_path)/focont"
	cd "$(sphinx_path)" && make markdown
	cp "$(sphinx_path)/build/markdown/focont.md" "$(doc_path)"

.PHONY: doc
