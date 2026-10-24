"""Templates of the indocyanine green models: the creators and the terms of use."""

from pkpdsim.templates import terms_of_use
from sbmlutils.factory import Creator

__all__ = ["creators", "terms_of_use"]

creators = [
    Creator(
        familyName="Köller",
        givenName="Adrian",
        email="adriankl39@googlemail.com",
        organization="Humboldt-University Berlin, Institute for Theoretical Biology",
        orcid="",
    ),
    Creator(
        familyName="Myshkina",
        givenName="Mariia",
        email="myshkina.mariia@gmail.com",
        organization="Humboldt-University Berlin, Institute for Theoretical Biology",
    ),
    Creator(
        familyName="König",
        givenName="Matthias",
        email="koenigmx@hu-berlin.de",
        organization="Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering",
        site="https://livermetabolism.com",
    ),
]
