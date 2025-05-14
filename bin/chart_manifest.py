#!/usr/bin/env python3

import json
import os

chart_manifest = []

for fp in os.listdir("."):
    if fp.endswith('.vt.json'):
        content = json.load(open(fp))
        chart_manifest.append(dict(
            type="vitessce",
            config=fp,
            name=content.get("name", ""),
            desc=content.get("description", "")
        ))

with open("chart.manifest.pipeline.json", "w") as handle:
    json.dump(chart_manifest, handle, indent=4)
