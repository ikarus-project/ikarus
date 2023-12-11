# SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

import os
import re
import datetime


def update_license_header(file_path):
    try:
        with open(file_path, "r", newline="") as file:
            content = file.read()

        # Update the license header using regex
        updated_content, nsubs = re.subn(
            r"SPDX-FileCopyrightText: (\d{4}-\d{4}) The Ikarus Developers",
            lambda match: f"SPDX-FileCopyrightText: 2021-{datetime.datetime.now().year} The Ikarus Developers",
            content,
        )
        if nsubs > 0:
            with open(file_path, "w", newline="") as file:
                file.write(updated_content)
            print(f"Updated license header in: {file_path}")
        else:
            print(f"Didn't update license header in: {file_path}")

    except Exception as e:
        print(f"Error updating license header in {file_path}: {e}")


def process_directory(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            update_license_header(file_path)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath("__file__"))

    print(script_dir)
    process_directory(script_dir)
