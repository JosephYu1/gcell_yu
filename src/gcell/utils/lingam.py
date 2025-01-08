# from cdt https://github.com/FenTechSolutions/CausalDiscoveryToolbox
import fileinput
import logging
import subprocess
import uuid
import warnings
from pathlib import Path
from shutil import copy, rmtree
from tempfile import gettempdir

import networkx as nx
from pandas import read_csv

RPATH = "Rscript"


def message_warning(msg, *a, **kwargs):
    """Ignore everything except the message."""
    return str(msg) + "\n"


warnings.formatwarning = message_warning
init = True


def launch_R_script(
    template, arguments, output_function=None, verbose=True, debug=False
):
    """Launch an R script, starting from a template and replacing text in file
    before execution."""
    base_dir = Path(gettempdir()) / f"cdt_R_script_{uuid.uuid4()!s}"
    base_dir.mkdir(parents=True, exist_ok=True)
    rpath = RPATH
    scriptpath = base_dir / f"instance_{Path(template).name}"
    copy(str(template), str(scriptpath))

    # Converting Paths to OS-compliant paths
    for arg in arguments:
        if isinstance(arguments[arg], Path | str):
            arguments[arg] = str(arguments[arg]).replace("\\", "\\\\")

    with fileinput.FileInput(str(scriptpath), inplace=True) as file:
        for line in file:
            mline = line
            for elt in arguments:
                mline = mline.replace(elt, arguments[elt])
            print(mline, end="")
    print(
        rpath,
        scriptpath,
        "Please make sure Rscript is in your PATH and pcalg is installed in R. Consider using mamba or conda to install r-pcalg.",
    )
    if output_function is None:
        try:
            output = subprocess.call(
                [str(rpath), "--no-restore --no-save --no-site-file", str(scriptpath)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except Exception as e:
            print("R Call errored, is R available ?")
            raise e

    else:
        try:
            if verbose:
                process = subprocess.Popen(
                    [
                        str(rpath),
                        "--no-restore --no-save --no-site-file",
                        str(scriptpath),
                    ]
                )
            else:
                process = subprocess.Popen(
                    [
                        str(rpath),
                        "--no-restore --no-save --no-site-file",
                        str(scriptpath),
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            process.wait()
        except KeyboardInterrupt:
            if not debug:
                rmtree(base_dir)
            raise KeyboardInterrupt
        try:
            output = output_function()
        # Cleaning up
        except Exception as e:
            if not debug:
                rmtree(base_dir)
            if not verbose:
                out, err = process.communicate()
                print("\nR Python Error Output \n-----------------------\n")
                print(e)
                raise RuntimeError(
                    "RProcessError \nR Process Error Output \n-----------------------\n"
                    + str(err, "ISO-8859-1")
                ) from None
            print("\nR Python Error Output \n-----------------------\n")
            print(e)
            raise RuntimeError("RProcessError ") from None

    if not debug:
        rmtree(base_dir)
    return output


class GraphModel:
    """Base class for all graph causal inference models.

    Usage for undirected/directed graphs and raw data. All causal discovery
    models out of observational data base themselves on this class. Its main
    feature is the predict function that executes a function according to the
    given arguments.
    """

    def __init__(self):
        """Init."""
        super().__init__()

    def predict(self, df_data, graph=None, **kwargs):
        """Orient a graph using the method defined by the arguments.

        Depending on the type of `graph`, this function process to execute
        different functions:

        1. If ``graph`` is a ``networkx.DiGraph``, then ``self.orient_directed_graph`` is executed.
        2. If ``graph`` is a ``networkx.Graph``, then ``self.orient_undirected_graph`` is executed.
        3. If ``graph`` is a ``None``, then ``self.create_graph_from_data`` is executed.

        Args:
            df_data (pandas.DataFrame): DataFrame containing the observational data.
            graph (networkx.DiGraph or networkx.Graph or None): Prior knowledge on the causal graph.

        .. warning::
           Requirement : Name of the nodes in the graph must correspond to the
           name of the variables in df_data
        """
        if graph is None:
            return self.create_graph_from_data(df_data, **kwargs)
        elif isinstance(graph, nx.DiGraph):
            return self.orient_directed_graph(df_data, graph, **kwargs)
        elif isinstance(graph, nx.Graph):
            return self.orient_undirected_graph(df_data, graph, **kwargs)
        else:
            print("Unknown Graph type")
            raise ValueError

    def orient_undirected_graph(self, data, umg, **kwargs):
        """Orient an undirected graph.

        .. note::
           Not implemented: will be implemented by the model classes.
        """
        raise NotImplementedError

    def orient_directed_graph(self, data, dag, **kwargs):
        """Re/Orient an undirected graph.

        .. note::
           Not implemented: will be implemented by the model classes.
        """
        raise NotImplementedError

    def create_graph_from_data(self, data, **kwargs):
        """Infer a directed graph out of data.

        .. note::
           Not implemented: will be implemented by the model classes.
        """
        raise NotImplementedError


class LiNGAM(GraphModel):
    r"""LiNGAM algorithm **[R model]**.


    **Description:** Linear Non-Gaussian Acyclic model. LiNGAM handles linear
    structural equation models, where each variable is modeled as
    :math:`X_j = \sum_k \alpha_k P_a^{k}(X_j) + E_j,  j \in [1,d]`,
    with  :math:`P_a^{k}(X_j)` the :math:`k`-th parent of
    :math:`X_j` and :math:`\alpha_k` a real scalar.

    **Required R packages**: pcalg

    **Data Type:** Continuous

    **Assumptions:** The underlying causal model is supposed to be composed of
    linear mechanisms and non-gaussian data. Under those assumptions, it is
    shown that causal structure is fully identifiable (even inside the Markov
    equivalence class).

    Args:
        verbose (bool): Sets the verbosity of the algorithm. Defaults to
           `cdt.SETTINGS.verbose`

    .. note::
       Ref: S.  Shimizu,  P.O.  Hoyer,  A.  Hyvärinen,  A.  Kerminen  (2006)
       A  Linear  Non-Gaussian  Acyclic Model for Causal Discovery;
       Journal of Machine Learning Research 7, 2003–2030.

    .. warning::
       This implementation of LiNGAM does not support starting with a graph.

    Example:
        >>> import networkx as nx
        >>> from cdt.causality.graph import LiNGAM
        >>> from cdt.data import load_dataset
        >>> data, graph = load_dataset("sachs")
        >>> obj = LiNGAM()
        >>> output = obj.predict(data)
    """

    def __init__(self, verbose=False):
        """Init the model and its available arguments."""
        logging.info("R Package pcalg is needed for LiNGAM.")

        super().__init__()

        self.arguments = {
            "{FOLDER}": "/tmp/cdt_LiNGAM/",
            "{FILE}": "/data.csv",
            "{VERBOSE}": "FALSE",
            "{OUTPUT}": "/result.csv",
        }
        self.verbose = "FALSE"

    def orient_undirected_graph(self, data, graph):
        """Run LiNGAM on an undirected graph."""
        # Building setup w/ arguments.
        raise ValueError("LiNGAM cannot (yet) be ran with a skeleton/directed graph.")

    def orient_directed_graph(self, data, graph):
        """Run LiNGAM on a directed_graph."""
        raise ValueError("LiNGAM cannot (yet) be ran with a skeleton/directed graph.")

    def create_graph_from_data(self, data):
        """Run the LiNGAM algorithm.

        Args:
            data (pandas.DataFrame): DataFrame containing the data

        Returns:
            networkx.DiGraph: Prediction given by the LiNGAM algorithm.

        """
        # Building setup w/ arguments.
        self.arguments["{VERBOSE}"] = str(self.verbose).upper()
        results = self._run_LiNGAM(data, verbose=self.verbose)

        return nx.relabel_nodes(
            nx.DiGraph(results), {idx: i for idx, i in enumerate(data.columns)}
        )

    def _run_LiNGAM(self, data, fixedGaps=None, verbose=True):
        """Setting up and running LiNGAM with all arguments."""
        # Run LiNGAM
        run_dir = Path(gettempdir()) / f"cdt_lingam_{uuid.uuid4()!s}"
        run_dir.mkdir(parents=True, exist_ok=True)
        self.arguments["{FOLDER}"] = str(run_dir)

        def retrieve_result():
            return read_csv(run_dir / "result.csv", delimiter=",").values

        try:
            data.to_csv(run_dir / "data.csv", header=False, index=False)
            template_path = Path(__file__).parent / "lingam.r"
            lingam_result = launch_R_script(
                template_path,
                self.arguments,
                output_function=retrieve_result,
                verbose=verbose,
            )
        # Cleanup
        except Exception as e:
            rmtree(run_dir)
            raise e
        except KeyboardInterrupt:
            rmtree(run_dir)
            raise KeyboardInterrupt
        rmtree(run_dir)
        return lingam_result


# %%
