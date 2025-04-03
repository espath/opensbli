#!/usr/bin/env python3

try:
    import jinja2
except ImportError:
    print("Trying to Install required module: jinja2\n")
    os.system('python -m pip install jinja2')

from jinja2 import Environment, FileSystemLoader
import os
try:
    import commentjson
except ImportError:
    print("Trying to Install required module: commentjson\n")
    os.system('python -m pip install commentjson')
import commentjson as json
import argparse

parser = argparse.ArgumentParser(
    prog = 'Submit.py',
    description = 'Automatically generate submission for a specified supercomputer.',
    epilog = '')
parser.add_argument('filename', type=str, nargs='?', default="JobSubmission.json",help="JSON file for specifying submission options.")

args = parser.parse_args()

with open(args.filename) as JsonOptionFile:
    JobOptions=json.load(JsonOptionFile)
JsonOptionFile.close()

environment = Environment(loader=FileSystemLoader("SubmissionTemplate/"))
ScriptTemplate = environment.get_template(JobOptions["Machine"])

Script = ScriptTemplate.render(JobOptions)
with open(JobOptions["JobName"], mode="w", encoding="utf-8") as message:
    message.write(Script)
