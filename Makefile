.PHONY: build test release clean install help
JULIA ?= julia

build:
	@echo "Julia is interpreted - no build step required"

test: build
	@$(JULIA) src/numerical_solver.jl

release:
	@tar czf dist/shae-julia-numerical-solver-v0.1.0-src.tar.gz src/ scripts/ Makefile README.md DESIGN.md LICENSE docs/

install:
	@echo "Run with: julia src/numerical_solver.jl"

clean:
	@rm -rf dist/

help:
	@echo "Targets: build test release install clean help"
