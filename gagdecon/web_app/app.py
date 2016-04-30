
import os

from flask import (
    Flask, request, session, g, redirect, url_for,
    abort, render_template, flash, Markup, make_response, jsonify,
    Response)


from jinja2 import FileSystemLoader


def prepare_env(env):
    env.loader = FileSystemLoader(os.path.join(os.path.dirname(__file__), 'html'))

app = Flask(__name__)

prepare_env(app.jinja_env)

DATABASE = None
DEBUG = True
SECRETKEY = 'LG0264AN0H5T5T1L5RUGVLW2FVUCDEMP'
SERVER = None


@app.route("/")
def index():
    return render_template("index.templ")


@app.route("/run", methods=["POST"])
def run_search():
    mass_error_tolerance = float(request.form["mass_error_tolerance"])
    max_charge = -abs(int(request.form['max_charge']))
    return redirect("/")


if __name__ == '__main__':
    app.run(port=5000, use_reloader=True)
