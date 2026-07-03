# module with utility functions for Pannotator
import gzip
import json
import pickle
from pathlib import Path


def load_pickle(pickle_path: Path) -> dict:
    with open(pickle_path, "rb") as f:
        return pickle.load(f)


def dump_pickle(data: dict, pickle_path: Path) -> dict:
    with open(pickle_path, "wb") as f:
        return pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_json(json_path: Path) -> dict:
    with open(json_path, "r") as f:
        return json.load(f)


def dump_json(data: dict, outpath: Path, indent: int | str | None = None) -> None:
    with open(outpath, "w") as f:
        json.dump(data, f, separators=(",", ":"), indent=indent)


def load_gzed_json(json_path: Path) -> dict:
    with gzip.open(json_path, "rt", encoding="utf-8") as f:
        return json.load(f)


def dump_gzed_json(data: dict, outpath: Path) -> None:
    with gzip.open(outpath, "wt", encoding="utf-8") as f:
        json.dump(data, f, separators=(",", ":"))


_FORMATS = {
    (".json",): (load_json, dump_json),
    (".json", ".gz"): (load_gzed_json, dump_gzed_json),
    (".pkl",): (load_pickle, dump_pickle),
    (".pickle",): (load_pickle, dump_pickle),
}


def _format(path: Path) -> tuple[str, ...]:
    if path.suffixes[-2:] == [".json", ".gz"]:
        return (".json", ".gz")
    return (path.suffix,)


def load_pangenome(path: Path):
    return _FORMATS[_format(path)][0](path)


def dump_pangenome(data: dict, path: Path):
    _FORMATS[_format(path)][1](data, path)
