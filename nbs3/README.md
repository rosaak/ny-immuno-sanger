# How to Run the script

- this script matches the heavey chain or light chain probes to abi files and matched files will be moved to a new dir.
- it also calcualtes the mean quality scores
- requires `fns.py`, `utils.py` and `config_single_vx.ini`

---

```bash
$ python3 vx-splitter-single.py --help
Usage: vx-splitter-single.py [OPTIONS] COMMAND [ARGS]...

Options:
  --install-completion  Install completion for the current shell.
  --show-completion     Show completion for the current shell, to copy it or
                        customize the installation.
  --help                Show this message and exit.

Commands:
  message
  single-dir-analysis
```

---

```bash
$ python3 vx-splitter-single.py message       
This script matches the abi file from a single dir with the probe file and move it to a different dir
```
---

```bash
$ python3 vx-splitter-single.py single-dir-analysis --help
Usage: vx-splitter-single.py single-dir-analysis [OPTIONS] CONFIG_FILE_PATH

Arguments:
  CONFIG_FILE_PATH  [required]

Options:
  --help  Show this message and exit.
```
  
---