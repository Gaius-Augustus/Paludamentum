# From Tiberius (52b87dd), tiberius/scripts/utility.py.
# Copyright (c) 2023 Lars Gabriel. MIT License, see LICENSE-Tiberius.
# Changed in Paludamentum: reduced to run_subprocess, the only function that
# hc_module.py uses.
import logging, subprocess

#Authors: "Amrei Knuth", "Lars Gabriel"
#Credits: "Katharina Hoff"
#Email:"lars.gabriel@uni-greifswald.de"
#Date: "Janurary 2025"

logger = logging.getLogger(__name__)


def run_subprocess(command, capture_output=True, text=True,
            check=True, error_message="Subprocess execution failed",
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, stdin=None):
    """
    Execute a subprocess command with error handling and logging.

    This function executes a given command using `subprocess.run`, captures its output,
    and handles errors gracefully. If the command fails, it logs the error and raises an exception.

    :param command: List of command arguments to execute (e.g., ["ls", "-l"]).
    :param capture_output: Whether to capture the command's stdout and stderr (default: True).
    :param text: If True, returns stdout and stderr as strings instead of bytes (default: True).
    :param check: If True, raises a `subprocess.CalledProcessError` on a non-zero return code (default: True).
    :param error_message: Custom error message to log and raise if the command fails (default: "Subprocess execution failed").
    :return: The result of the subprocess (subprocess.CompletedProcess object).
    :raises RuntimeError: If the subprocess fails and `check` is set to False.
    :raises subprocess.CalledProcessError: If the subprocess fails and `check` is set to True.
    """
    try:
        logging.info(f"Executing command: {' '.join(command)}")
        if capture_output:
            result = subprocess.run(command, capture_output=capture_output, text=text,
                            check=check, shell=shell, stdin=stdin)
        else:
            result = subprocess.run(command, capture_output=capture_output, text=text, check=check,
                            stdout=stdout, stderr=stderr, shell=shell, stdin=stdin)
        if result.returncode == 0:
            logging.info("Subprocess completed successfully.")
        else:
            logging.warning(f"Subprocess completed with non-zero return code: {result.returncode}")

        return result

    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message}: {e}")
        logging.error(f"Command: {' '.join(command)}")
        logging.error(f"stdout: {e.stdout}")
        logging.error(f"stderr: {e.stderr}")
        raise

    except Exception as e:
        logging.error(f"Unexpected error while running subprocess: {e}")
        logging.error(f"Command: {' '.join(command)}")
        raise RuntimeError(f"{error_message}: {e}") from e


