# -*- coding: utf-8 -*-
from sqlalchemy_utils import database_exists, create_database
from project_spectra.models import engine
import os
from pathlib import Path

# Init Database
# If it doesn't exist, create one
if not database_exists(engine.url):
    create_database(engine.url)

# identify home directory
home_dir = Path.home()

# join paths
data_path = os.path.join(str(home_dir), ".group_project", "danqi", "../../data")
result_path = os.path.join(str(home_dir), ".group_project", "danqi", "result")

# create folder
os.makedirs(data_path, mode=0o777, exist_ok=True)
os.makedirs(result_path, mode=0o777, exist_ok=True)
