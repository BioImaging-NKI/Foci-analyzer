from .cli import main

if __name__ == "__main__":  # the batch run starts worker processes that import this file again
    raise SystemExit(main())
