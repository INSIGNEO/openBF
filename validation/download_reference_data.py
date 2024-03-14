import requests
from pathlib import Path


def download_reference_data(simulation, gist_id):
    if not (p := Path(simulation, f"{simulation}_ref")).exists():
        p.mkdir()

    response = requests.get(f"https://api.github.com/gists/{gist_id}")

    if response.status_code == 200:
        gist_data = response.json()

        for filename, file_info in gist_data["files"].items():
            if (ref_file := Path(p, filename)).exists():
                print(f"{ref_file} already there, skipping download")
                continue

            file_content = requests.get(file_info["raw_url"]).text

            with ref_file.open("w") as f:
                f.write(file_content)

            print(f"{ref_file} downloaded")
    else:
        raise FileNotFoundError(
            f"Failed to retrieve Gist: {response.status_code} - {response.text}"
        )


if __name__ == "__main__":
    folders = ["cca", "uta", "ibif"]
    gist_ids = [
        "d0075e8d50ef6eefc800a1fa4545b2b1",
        "6cff8b76b77d2c405a16fa9e990d14a7",
        "66dfe0f8b0343982c9e0a9edaef746f4",
    ]

    for folder, gist_id in zip(folders, gist_ids):
        download_reference_data(folder, gist_id)
