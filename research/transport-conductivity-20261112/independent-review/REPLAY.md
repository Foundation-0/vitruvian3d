Run the independent arithmetic checker from the transport note directory:

```sh
python3 independent-review/independent_bounds.py --source-root /path/to/vitruvian3d
```

The source checkout must contain the four files pinned to commit `e7445d56b9ecfd97b408653a2726f3bab70a6e4c`. The checker verifies their working-tree bytes and their bytes at that Git commit. It also checks the frozen theorem, author checker, and certificate beside this directory.

The default output is `independent-review/independent-bounds.json`. Use `--output /path/to/result.json` to write elsewhere. The report includes the checkout's current head, which may be later than the pinned source commit when the four required files are unchanged.
