"""
Description
-----------

Parses templatized GitHub Issue into a DataClass and dumps to JSON.

"""

import json
import os
import sys
from dataclasses import dataclass, asdict

@dataclass
class NewProfilingSite:
    """Dataclass representing contents of templatized GitHub Issue

    Attributes:
        name (str): Name given to new profiling site
        latitude (float): Latitude for new profiling site (deg N)
        longitude (float): Longitude for new profiling site (deg E)
        models (list): List of models new profiling site will be added to
    Methods:
        None

    Usage:
        # Create an instance of NewProfilingSite
        new_profiling_site = NewProfilingSite(issue_body)

        # Access the attributes of the instance
        print(new_profiling_site.<attribute>)
    """

    name: str
    latitude: float
    longitude: float
    models: list[str]

    def __init__(self, issue_body: str) -> None:
        """Initialize the NewProfilingSite instance

        Performs tasks such as parsing GitHub Issue body, identifying the entries, identifying
        their values, creating a dict from these entries/values, and setting attributes of
        dataclass.

        :param issue_body: templatized Github Issue body
        :type issue_body: str
        :returns: {None}
        :rtype: {None}
        """

        # split raw string by new lines, drop empty strings
        keys_and_values = [entry for entry in issue_body.split("\n") if entry != ""]

        # keys lead with "### "
        keys = [key.replace("### ", "") for key in keys_and_values if "###" in key]

        # if no leading "### " and not blank, then value
        values = [value for value in keys_and_values if "###" not in value]

        # convert keys and values to dict
        issue_dict = {key.replace(" ", ""): value for key,value in zip(keys,values)}

        # required fields:
        required_fields = [
            "ProfilingSiteName",
            "ProfilingSiteLatitude",
            "ProfilingSiteLongitude",
            "Models"]

        # verify that we have all required fields in issue_dict before proceeding. If not, exit.
        if not all(field in issue_dict for field in required_fields):
            sys.exit()

        # set class attributes
        self.name = issue_dict["ProfilingSiteName"]
        self.latitude = float(issue_dict["ProfilingSiteLatitude"])
        self.longitude = float(issue_dict["ProfilingSiteLongitude"])
        self.models = [model.lower() for model in issue_dict["Models"].split(", ")]

def to_json(data: NewProfilingSite, file_name: str) -> None:
    """Dump contents of dataclass to JSON

    :param data: dataclass containing contents of templatized GitHub issue
    :type data: dataclass
    :param file_name: JSON file name to write templatized GitHub issue to
    :type file_name: str
    :returns: {None}
    :rtype: {None}
    """

    with open(file_name, "w", encoding="utf-8") as issue_json:
        json.dump(asdict(data), issue_json)

def main() -> None:
    """Main method responsible for executing all steps of script

    Script flows as follows:

    #. If templatized GitHub Issue exists as environment variable "body",
       retrieve it.
    #. Parse issue body into dataclass.
    #. Write parsed issue body to JSON from dataclass.

    :param {None}: {None}
    :returns: {None}
    :rtype: {None}
    """

    # get issue body via BODY environment variable
    if issue_body := os.environ.get("body"):

        # instantiate NewProfilingSite dataclass from issue body
        new_profiling_site = NewProfilingSite(issue_body)

        # dump profiling site characteristics to json
        to_json(new_profiling_site, "new_profiling_site.json")

# run main routine
if __name__ == '__main__':
    main()
