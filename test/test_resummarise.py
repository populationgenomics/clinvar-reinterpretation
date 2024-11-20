from copy import deepcopy
from datetime import datetime

import pytest
import zoneinfo

from clinvarbitration.resummarise_clinvar import Consequence, Submission, consequence_decision

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
BASIC_SUB = Submission(datetime.now(tz=TIMEZONE), 'foo', Consequence.UNKNOWN, 'review')


@pytest.mark.parametrize(
    ('consequences', 'expected'),
    [
        ([Consequence.PATHOGENIC], Consequence.PATHOGENIC),
        ([Consequence.PATHOGENIC, Consequence.PATHOGENIC], Consequence.PATHOGENIC),
        ([Consequence.PATHOGENIC, Consequence.UNCERTAIN], Consequence.PATHOGENIC),
        ([Consequence.BENIGN], Consequence.BENIGN),
        ([Consequence.BENIGN, Consequence.PATHOGENIC], Consequence.CONFLICTING),
    ],
)
def test_consequence_decision_path_single(
    consequences: list[Consequence],
    expected: Consequence,
):
    all_subs = []
    for con in consequences:
        sub = deepcopy(BASIC_SUB)
        sub.classification = con
        all_subs.append(sub)
    assert consequence_decision(all_subs, False) == expected
