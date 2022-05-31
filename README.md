
### clone this repo

```bash
git clone https://github.com/rosaak/ny-immuno-sanger.git
cd ny-immuno-sanger 
```

### setting up the env

```bash
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

### file structure

- Excel File
```raw
14 clones H3.xlsx

```

- AB1 files
```raw
Sequencing results
├── VH
│   ├── FST470-hIgG1-011-A10_A10-VH60.ab1
│   ├── FST470-hIgG1-011-A11_A11-VH60.ab1
│   ├── FST470-hIgG1-011-A12_A12-VH60.ab1
│   ├── FST470-hIgG1-011-A1_A01-VH60.ab1
...
└── VL
    ├── FST470-hIgG1-011-A10_A10-VL79.ab1
    ├── FST470-hIgG1-011-A11_A11-VL79.ab1
    ├── FST470-hIgG1-011-A12_A12-VL79.ab1
...

```

- Template
```raw
Templates
├── VH
│   ├── VH-FST470-hF-010-A01.gb
│   ├── VH-FST470-hF-010-A02.gb
│   ├── VH-FST470-hF-010-A05.gb
...
└── VL
    ├── VL-FST470-hF-010-A01.gb
    ├── VL-FST470-hF-010-A02.gb
    ├── VL-FST470-hF-010-A05.gb

```

### what is happening
- pick the pattern_sequence from the Excel file and search it against files in `Sequencing results/{VH,VL}/*.ab1`.
- if there is a match then create a new dir `Match_Res/{pattern_name}` and put the matched `.ab1` in it into VH and VL  

