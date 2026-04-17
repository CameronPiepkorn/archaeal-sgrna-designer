# Contributing to archaeal-sgrna-designer

Thank you for your interest in contributing! This document outlines how to get started.

## Reporting bugs

Please open an issue on GitHub with:
- A description of the bug
- A minimal reproducible example
- Your Python version and operating system

## Suggesting features

Open an issue describing:
- The feature you would like
- The biological use case it addresses
- Any relevant literature

## Contributing code

1. Fork the repository
2. Create a new branch: `git checkout -b my-feature`
3. Make your changes
4. Run the tests: `pytest`
5. Push your branch and open a pull request

## Development setup

```bash
git clone https://github.com/CameronPiepkorn/archaeal-sgrna-designer
cd archaeal-sgrna-designer
pip install -e ".[dev]"
pytest
```

## Code style

- Follow PEP 8
- Add docstrings to new functions
- Add tests for new functionality
- Keep biological claims grounded in published literature with citations

## Questions

Open an issue or start a discussion on GitHub.