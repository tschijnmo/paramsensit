"""
Test the sensitivity of simulation parameters

This is a very small package for testing the sensitivity of simulation
results on some parameters. Its main functions include

.. autosummary::
    :toctree:

    run_simuls
    get_results

"""


import os
import subprocess
import collections

from jinja2 import Environment, FileSystemLoader, TemplateNotFound


#
# Public functions
# ----------------
#


def run_simuls(params, files, acts):
    """Performs the simulations to tests the sensitivity

    This is the driver function for running the simulations before gathering
    the results.

    :param params: The sequence of parameters whose sensitivity is to be tested.
        The values need to be a triple of parameter name, parameter value,
        and a differentiation step size for changing the parameter.
    :param files: An iterable of the names of the input files of the simulation.
        All the files should be present in the current working directory and
        are going to be written into the subdirectories for the simulations.
    :param acts: An iterable of actions that are going to be performed in each
        of the subdirectories to initiate the simulations. Needs to be given as
        strings for the command that is to be run.
    :returns: 0 for success.
    """

    # Set the Jinja template engine up and get all the templates.
    orig_cwd = os.getcwd()
    env = Environment(loader=FileSystemLoader(orig_cwd))
    templs = []
    for i in files:
        try:
            templs.append(env.get_template(i))
        except TemplateNotFound:
            raise ValueError(
                'Input file template {} cannot be found'.format(i)
            )
        continue

    # Do the job for each of the given parameters.
    for idx, param in enumerate(params):

        # First we need to get all the parameter values and sub-directory names.
        param_vals = _get_param_vals(param)
        dir_names = _get_dir_names(param)

        # Render the templates and run the requested actions in each of the
        # sub-directories.

        for param_val, dir_name in zip(param_vals, dir_names):

            # Make and switch to the subdirectory.
            os.makedirs(dir_name, exist_ok=True)
            os.chdir(os.path.join(orig_cwd, dir_name))
            print(
                'Parameter {p_name} value {val} in directory {dir_name}'.format(
                    p_name=param[0], val=param_val, dir_name=dir_name
                )
            )

            # Make the template rendering context.
            ctx = {name: val for name, val, _ in params}
            ctx[param[0]] = param_val
            ctx['dir_name'] = dir_name

            # Generate all the templates.
            for templ, name in zip(templs, files):
                with open(name, 'w') as input_fp:
                    input_fp.write(templ.render(ctx))
                continue

            # Make all the requested actions.
            for i in acts:
                subprocess.Popen(i.split())
                continue

            # Switch back to the original working directory.
            os.chdir(orig_cwd)

            # Continue to the next value of the parameter
            continue

        # Continue to the next parameter.
        continue

    return 0


ParamSensit = collections.namedtuple(
    'ParamSensit', [
        'value',
        'deriv',
        'cond',
    ]
)


def get_results(params, acts):
    """Gets the results about the sensitivity

    :param params: An iterable of the parameters whose sensitivity are to be
        tested, given as triples of name, value, and difference.
    :param acts: The actions to be performed for the results. It should be a
        mapping with names of the results as keys and callables for getting
        the results as values. The callables are going to be called with no
        arguments in the subdirectories for the simulations.
    :returns: A dictionary giving the results for each of the parameters on
        each of the results. With the pair of the name of the parameter and
        the name of the result as keys and ``ParamSensit`` objects as values.
    :rtype: dict
    """

    orig_cwd = os.getcwd()
    res = {}

    # Make the processing for each of the parameters one-by-one.
    for param in params:

        # The directory names and the values.
        dir_names = _get_dir_names(param)

        # Get the simulation results one-by-one in the order.
        res_vals = {i: [] for i in acts.keys()}
        for dir_name in dir_names:

            # Switch the directory.
            os.chdir(os.path.join(orig_cwd, dir_name))

            for k, v in acts.items():
                res_vals[k].append(
                    v()
                )
                continue

            # Continue to the next parameter value/subdirectory.
            os.chdir(orig_cwd)
            continue

        # Compute the results for the current parameter.
        res.update(
            _compute_results(param, res_vals)
        )

        # Continue to the next parameter.
        continue

    # Print the results on the stdout before returning.
    _print_res(res)
    return res


#
# Private functions
# -----------------
#
# Directory and values generation
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#


def _get_param_vals(param):
    """Gets the value of the parameters

    :param param: The parameter triple whose values are to be generated.
    :returns: The value of the parameter on a five-point stencil.
    """

    return [
        param[1] + i * param[2]
        for i in range(-2, 3)
    ]


def _get_dir_names(param):
    """Gets the names of the directories for testing the given parameter

    :param param: The parameter triple to be tested.
    :returns: A list of directories names to be used to test the given
        parameter.
    """

    return [
        '-'.join([param[0], str(i)])
        for i in range(0, 5)
    ]


#
# Result computation
# ^^^^^^^^^^^^^^^^^^
#


def _compute_results(param, res_vals):
    """Computes the results for the sensitivity of a parameter

    :param param: The parameter whose sensitivity is to be tested.
    :param res_vals: The dictionary giving the results of the simulation at
        different values of the parameter. The names of the results are keys
        and the five values of the results are the values.
    :returns: A dictionary giving the sensitivity of the results to the given
        parameter.
    """

    res = {}

    # Treat the results one-by-one.
    for name, vals in res_vals.items():
        res[(param[0], name)] = _compute_result(param, vals)
        continue

    return res


def _compute_result(param, vals):
    """Computes the sensitivity result for a single parameter on a single result

    :param param: The parameter to be tested.
    :param vals: A list of five points for the results.
    :returns: A parameter sensitivity for this pair.
    """

    # The simulation result at the given parameter value.
    value = vals[2]

    # The five-point stencil derivative formula.
    deriv = (
        -vals[4] + 8 * vals[3] - 8 * vals[1] + vals[0]
    ) / (12 * param[2])

    # The condition number.
    cond = (param[1] * deriv) / value

    return ParamSensit(
        value=value, deriv=deriv, cond=cond
    )


#
# Result printing
# ^^^^^^^^^^^^^^^
#


def _print_res(res):
    """Prints the results on the screen

    :param res: The dictionary of the results.
    :returns: 0 for success.
    """

    tot_length = 120
    n_fields = 5
    fmt = ''.join(['{:^', str(tot_length // n_fields), '}']) * n_fields

    print(fmt.format(
        'Parameter', 'Result', 'Value', 'Derivative', 'Condition N'
    ))
    for k, v in res.items():
        print(
            fmt.format(*(k + tuple(v)))
        )
        continue

    return 0
