#!/bin/bash

echo "Installing or updating Ladybug/Honeybee"

pip install lbt-ladybug --upgrade
pip install ladybug-comfort --upgrade
pip install ladybug-geometry --upgrade

pip install lbt-dragonfly --upgrade
pip install uwg --upgrade

pip install lbt-honeybee --upgrade
pip install honeybee-core --upgrade
pip install honeybee-schema --upgrade
pip install honeybee-radiance --upgrade
pip install honeybee-energy --upgrade
pip install honeybee-energy-standards --upgrade