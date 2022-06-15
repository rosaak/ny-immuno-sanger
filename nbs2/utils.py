import logging
from datetime import datetime as dt
from pathlib import Path


def now():
    now = dt.now()
    return(now.strftime("%Y%m%d%H%M%S"))

def make_logger(name: str):
    """
    Make a logger
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create a file handler
    handler = logging.FileHandler(f"{name}_{now()}.log")
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return(logger)

def dir_exists_and_create(dir_path:str, create_dir:False) -> int:
    """
    return status
    1: dir_path created
    2: dir_path not exist and didn't create
    3: dir_path exist and didn't create
    4: dir_path already exists
    """
    if not Path(dir_path).exists() and create_dir:
        Path(dir_path).resolve().mkdir(parents=True, exist_ok=True)
        return(1)
    elif not Path(dir_path).exists() and not create_dir:
        return(2)
    elif Path(dir_path).exists() and not create_dir:
        return(3)
    elif Path(dir_path).exists() and create_dir:
        return(4)
    
def file_exists(file_path:str) -> int:
    """
    returns 
    1 : file exists
    0 : file doesn't exists
    """
    if Path(file_path).exists():
        return 1
    else:
        return 0

def main():
    pass

if __name__ == '__main__':
    main()