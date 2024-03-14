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


GISTS = {
    "cca": "11645e64d3583fee37f1b48e5a897b51",
    "uta": "6cff8b76b77d2c405a16fa9e990d14a7",
    "ibif": "feaae8d7d29cc9fa179e109f68abc08",
}
if __name__ == "__main__":
    for simulation, gist_id in GISTS.items():
        download_reference_data(simulation, gist_id)
