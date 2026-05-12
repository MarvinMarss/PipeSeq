# Release Checklist

Use this checklist before pushing a public release.

1. Confirm `Script/settings.json` contains no private paths before sharing screenshots or logs.
2. Confirm generated logs and results are ignored.
3. Run:

```bash
python -m py_compile Script/*.py
```

4. Launch the GUI with:

```bash
python Script/PipeSeq.py
```

5. Verify README installation steps against a clean environment.
6. Tag the release after pushing:

```bash
git tag vX.Y.Z
git push origin vX.Y.Z
```
