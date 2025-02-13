from copy import deepcopy
from datetime import datetime

import pytest
import zoneinfo

from clinvarbitration.scripts.resummarise_clinvar import Consequence, Submission, check_stars, consequence_decision

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
BASIC_SUB = Submission(datetime.now(tz=TIMEZONE), 'foo', Consequence.UNKNOWN, 'review')


def test_check_stars_none():
    """
    check that we always get the right number of stars!
    """
    # Unknown, should be skipped over
    dud_1 = deepcopy(BASIC_SUB)
    dud_2 = deepcopy(BASIC_SUB)
    # Pathogenic, but has a no-star rating
    dud_2.classification = Consequence.PATHOGENIC
    dud_2.review_status = 'no assertion criteria provided'
    assert check_stars([dud_1, dud_2]) == 0


def test_check_stars_1():
    """
    check that we always get the right number of stars!
    """
    expect_1_path = deepcopy(BASIC_SUB)
    expect_1_path.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path]) == 1
    expect_1_benign = deepcopy(BASIC_SUB)
    expect_1_benign.classification = Consequence.BENIGN
    assert check_stars([expect_1_benign]) == 1
    assert check_stars([expect_1_path, expect_1_benign]) == 1


def test_check_stars_3():
    """
    check that we always get the right number of stars!
    """
    expect_1_path = deepcopy(BASIC_SUB)
    expect_1_path.classification = Consequence.PATHOGENIC
    expect_1_neutral = deepcopy(BASIC_SUB)
    assert check_stars([expect_1_path, expect_1_neutral]) == 1
    expert_panel = deepcopy(BASIC_SUB)
    expert_panel.review_status = 'reviewed by expert panel'
    expert_panel.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path, expect_1_neutral, expert_panel]) == 3  # noqa: PLR2004


def test_check_stars_4():
    """
    check that we always get the right number of stars!
    """
    expect_1_path = deepcopy(BASIC_SUB)
    expect_1_path.classification = Consequence.PATHOGENIC
    expect_1_neutral = deepcopy(BASIC_SUB)
    assert check_stars([expect_1_path, expect_1_neutral]) == 1
    expert_panel = deepcopy(BASIC_SUB)
    expert_panel.review_status = 'reviewed by expert panel'
    expert_panel.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path, expect_1_neutral, expert_panel]) == 3  # noqa: PLR2004
    practice_guideline = deepcopy(BASIC_SUB)
    practice_guideline.review_status = 'practice guideline'
    practice_guideline.classification = Consequence.PATHOGENIC
    assert check_stars([expect_1_path, expect_1_neutral, expert_panel, practice_guideline]) == 4  # noqa: PLR2004


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
    assert consequence_decision(all_subs) == expected
