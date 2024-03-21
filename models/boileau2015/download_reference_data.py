import sys
import tempfile
import time
import zipfile
from pathlib import Path

import requests

SIMULATIONS = ["cca", "uta", "ibif", "adan56"]


def make_folder_structure():
    for simulation in SIMULATIONS:
        if not (p := Path(simulation, f"{simulation}_ref")).exists():
            p.mkdir()


def download_boileau2015_data(dst_folder):
    retry = 5
    while retry > 1:
        print("Downloading Boileau2015 supplementary material")
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
        }
        response = requests.get(
            "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fcnm.2732&file=cnm2732-sup-0001-Supplementary.zip",
            headers=headers,
        )

        if response.status_code == 200:
            with tempfile.TemporaryDirectory() as temp_dir:
                zip_path = Path(temp_dir, "nm2732-sup-0001-Supplementary.zip")
                with zip_path.open("wb") as z:
                    z.write(response.content)

                print("Extracting files")
                with zipfile.ZipFile(zip_path, "r") as z:
                    z.extractall(temp_dir)

                print("Copying files")
                root_path = Path(temp_dir, "SuppMaterial_1DBenchmark_Boileau_etal")
                for sim, folder in zip(
                    SIMULATIONS,
                    [
                        "Benchmark2_CommonCarotidArtery",
                        "Benchmark3_UpperThoracicAorta",
                        "Benchmark4_AorticBifurcation",
                    ],
                ):
                    print(sim)
                    tmp_data = Path(root_path, folder, "NumData")
                    dst_path = Path(dst_folder, sim, f"{sim}_ref")

                    if len(list(dst_path.glob("*"))) == 0:
                        tmp_data.rename(dst_path)
            sys.exit()

        else:
            print("Failed to retrieve Boileau2015 data, retrying...")
            retry -= 1
            time.sleep(5)

    print("Failed downloading data after 5 retries")


if __name__ == "__main__":
    make_folder_structure()
    download_boileau2015_data(Path.cwd())
