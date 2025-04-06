# 1) SETUP & IMPORTS

from ortools.sat.python import cp_model
import pandas as pd
from collections import defaultdict
from flask import Flask, request, jsonify
# from pyngrok import ngrok
from flask_cors import CORS
# -------------------------
# Global default values
days = ["Mon", "Tue", "Wed", "Thu", "Fri"]
timeslots = [
    "08:30-9:20",
    "09:20-10:10",
    "10:10-11:00",
    "11:30-12:20",
    "12:20-1:10",
    "2:00-2:50",
    "2:50-3:40",
    "3:40-4:30",
]
num_timeslots = len(timeslots)

# Globals to be set from API input.
room_info = {}         # { room_id: {"floor": ..., "has_speaker": ..., "has_multimedia": ...} }
lab_room_info = {}     # similarly for lab rooms
theory_rooms_original = []  # list(room_info.keys())
lab_rooms = []         # list(lab_room_info.keys())

# 2) HELPER FUNCTIONS
def subcourse_length(sc_code):
    if "_T1Sep" in sc_code:
        return 1
    elif "_T2Cons" in sc_code or "_L2Cons" in sc_code:
        return 2
    elif "_T4" in sc_code:
        return 4
    return 0

def expand_course(course):
    out = []
    nm = course["name"]
    ch = course["credit_hours"]
    is_lb = course["is_lab"]
    cid = course.get("id")
    if not is_lb:
        if ch == 2:
            out.append({
                "code": f"{nm}_T2Cons",
                "theory_hours": 2,
                "needs_lab_2slots": True,
                "is_lab_sub": False,
                "course_id": cid
            })
        elif ch == 3:
            out.append({
                "code": f"{nm}_T2Cons",
                "theory_hours": 2,
                "needs_lab_2slots": True,
                "is_lab_sub": False,
                "course_id": cid
            })
            out.append({
                "code": f"{nm}_T1Sep",
                "theory_hours": 1,
                "needs_lab_2slots": False,
                "is_lab_sub": False,
                "course_id": cid
            })
        elif ch == 4:
            out.append({
                "code": f"{nm}_T4",
                "theory_hours": 4,
                "needs_lab_2slots": False,
                "is_lab_sub": False,
                "course_id": cid
            })
    else:
        if ch == 3:
            out.append({
                "code": f"{nm}_T2Cons",
                "theory_hours": 2,
                "needs_lab_2slots": True,
                "is_lab_sub": False,
                "course_id": cid
            })
            out.append({
                "code": f"{nm}_L2Cons",
                "theory_hours": 0,
                "needs_lab_2slots": True,
                "is_lab_sub": True,
                "course_id": cid
            })
        elif ch == 4:
            out.append({
                "code": f"{nm}_T2Cons",
                "theory_hours": 2,
                "needs_lab_2slots": True,
                "is_lab_sub": False,
                "course_id": cid
            })
            out.append({
                "code": f"{nm}_T1Sep",
                "theory_hours": 1,
                "needs_lab_2slots": False,
                "is_lab_sub": False,
                "course_id": cid
            })
            out.append({
                "code": f"{nm}_L2Cons",
                "theory_hours": 0,
                "needs_lab_2slots": True,
                "is_lab_sub": True,
                "course_id": cid
            })
    return out

def is_first_third_year(y):
    return (y == 1 or y == 3)

def is_second_fourth_year(y):
    return (y == 2 or y == 4)

# 3) MODEL BUILDING FUNCTION
def build_model_no_gap_with_teachers(
    disc_info,
    teacher_info,
    locked_slots=None,     # list of dicts
    disabled_days=None,    # {(disc, year, section_index): [day, ...], ...}
    skip_day_pref_for=None,
    skip_single_room_for=None
):
    if locked_slots is None:
        locked_slots = []
    if disabled_days is None:
        disabled_days = {}
    if skip_day_pref_for is None:
        skip_day_pref_for = set()
    if skip_single_room_for is None:
        skip_single_room_for = set()

    model = cp_model.CpModel()
    # all_sections stores tuples: (disc, year, section_index, batch_id, section_id)
    all_sections = []
    theory_room_assign = {}
    lab_room_assign = {}
    schedule = {}
    course_subs = {}
    theory_usage = {}
    lab_usage = {}
    teacher_assign = {}

    # Normalize teacher time preferences.
    teacher_time_pref_list = {}
    for tID, tdata in teacher_info.items():
        prefs = tdata.get("time_prefs", [])
        if isinstance(prefs, dict):
            prefs = [prefs]
        # Convert "batch" to a tuple if it's a list.
        for pref in prefs:
            if "batch" in pref and isinstance(pref["batch"], list):
                pref["batch"] = tuple(pref["batch"])
        teacher_time_pref_list[tID] = prefs


    # Use year and section index for slot placement.
    for (disc, year), info in disc_info.items():
        sections_list = info["sections"]  # actual section IDs for response
        batch_id = info["batch_id"]       # preserved for response
        for idx, sec in enumerate(sections_list):
            section_index = idx
            all_sections.append((disc, year, section_index, batch_id, sec))
            for r in theory_rooms_original:
                theory_room_assign[(disc, year, section_index, r)] = model.NewBoolVar(f"ThRoom_{disc}_{year}_{section_index}_r{r}")
            expanded = []
            for c in info["courses"]:
                expanded.extend(expand_course(c))
            course_subs[(disc, year, section_index)] = expanded
            for sc in expanded:
                scode = sc["code"]
                for d in days:
                    for t_i in range(num_timeslots):
                        schedule[(disc, year, section_index, scode, d, t_i)] = model.NewBoolVar(f"Sched_{disc}_{year}_{section_index}_{scode}_{d}_{t_i}")
                if sc["is_lab_sub"]:
                    for lr in lab_rooms:
                        lab_room_assign[(disc, year, section_index, scode, lr)] = model.NewBoolVar(f"LabRm_{disc}_{year}_{section_index}_{scode}_r{lr}")
    for lock in locked_slots:
        d = lock['disc']
        yr = lock['year']    # now year (e.g., 1,2,…)
        sec_idx = lock['section']  # now section index (0,1,…)
        scode = lock['scode']
        day = lock['day']
        t_i = lock['timeslot']
        model.Add(schedule[(d, yr, sec_idx, scode, day, t_i)] == 1)
    for (d, yr, sec_idx, _, _ ) in [ (disc, year, section_index, batch_id, sec) for (disc, year, section_index, batch_id, sec) in all_sections]:
        disabled_day_list = disabled_days.get((d, yr, sec_idx), [])
        for day in disabled_day_list:
            for sc in course_subs[(d, yr, sec_idx)]:
                scode = sc["code"]
                for t_i in range(num_timeslots):
                    model.Add(schedule[(d, yr, sec_idx, scode, day, t_i)] == 0)

    valid_teachers_for = defaultdict(list)
    for tID, tdata in teacher_info.items():
        for assignment in tdata["assignments"]:
            # Expecting [discipline, batchId, year, sectionIndex, course_name, teacher_type]
            if len(assignment) == 6:
                di, batch_id, yr, sx, baseCourse, ttype = assignment
            elif len(assignment) == 5:
                # If only 5 values, assume a default year (e.g., 1)
                di, batch_id, sx, baseCourse, ttype = assignment
                yr = 1  # You can change this default if needed.
            else:
                raise ValueError(f"Unexpected teacher assignment format: {assignment}")
            if (di, yr, sx) not in course_subs:
                continue
            for sc in course_subs[(di, yr, sx)]:
                scode = sc["code"]
                baseName = scode.split("_")[0]
                sub_is_lab = sc["is_lab_sub"]
                if baseName == baseCourse:
                    if ttype == "lab" and sub_is_lab:
                        valid_teachers_for[(di, yr, sx, scode)].append(tID)
                    elif ttype == "theory" and not sub_is_lab:
                        valid_teachers_for[(di, yr, sx, scode)].append(tID)
                    elif ttype == "both":
                        valid_teachers_for[(di, yr, sx, scode)].append(tID)


    for (disc, yr, sec_idx, _, _) in all_sections:
        for sc in course_subs[(disc, yr, sec_idx)]:
            scode = sc["code"]
            tlist = valid_teachers_for.get((disc, yr, sec_idx, scode), [])
            if tlist:
                tvs = []
                for tID in tlist:
                    tid_str = str(tID)
                    tv = model.NewBoolVar(f"Tassign_{disc}_{yr}_{sec_idx}_{scode}_{tid_str}")
                    tvs.append(tv)
                    teacher_assign[(disc, yr, sec_idx, scode, tid_str)] = tv
                model.Add(sum(tvs) == 1)

    subUsedVars_for = defaultdict(list)
    for (di, yr, sec_idx, scode, tID) in teacher_assign:
        tv = teacher_assign[(di, yr, sec_idx, scode, tID)]
        for d in days:
            for t_i in range(num_timeslots):
                x = model.NewBoolVar(f"Used_{tID}_{di}_{yr}_{sec_idx}_{scode}_{d}_{t_i}")
                model.Add(x <= schedule[(di, yr, sec_idx, scode, d, t_i)])
                model.Add(x <= tv)
                model.Add(x >= schedule[(di, yr, sec_idx, scode, d, t_i)] + tv - 1)
                subUsedVars_for[(tID, d, t_i)].append(x)
    for (tID, d, t_i), usage_list in subUsedVars_for.items():
        model.Add(sum(usage_list) <= 1)

    for (disc, yr, sec_idx, _, _) in all_sections:
        exs = course_subs[(disc, yr, sec_idx)]
        for sc in exs:
            scode = sc["code"]
            if sc["is_lab_sub"]:
                for d in days:
                    for t_i in range(num_timeslots):
                        for lr in lab_rooms:
                            z = model.NewBoolVar(f"zLab_{disc}_{yr}_{sec_idx}_{scode}_{d}_{t_i}_{lr}")
                            model.Add(z <= schedule[(disc, yr, sec_idx, scode, d, t_i)])
                            model.Add(z <= lab_room_assign[(disc, yr, sec_idx, scode, lr)])
                            model.Add(z >= schedule[(disc, yr, sec_idx, scode, d, t_i)] + lab_room_assign[(disc, yr, sec_idx, scode, lr)] - 1)
                            lab_usage[(disc, yr, sec_idx, scode, d, t_i, lr)] = z
            else:
                for d in days:
                    for t_i in range(num_timeslots):
                        for r in theory_rooms_original:
                            z = model.NewBoolVar(f"zTheory_{disc}_{yr}_{sec_idx}_{scode}_{d}_{t_i}_{r}")
                            model.Add(z <= schedule[(disc, yr, sec_idx, scode, d, t_i)])
                            model.Add(z <= theory_room_assign[(disc, yr, sec_idx, r)])
                            model.Add(z >= schedule[(disc, yr, sec_idx, scode, d, t_i)] + theory_room_assign[(disc, yr, sec_idx, r)] - 1)
                            theory_usage[(disc, yr, sec_idx, scode, d, t_i, r)] = z

    for d in days:
        for t_i in range(num_timeslots):
            for r in theory_rooms_original:
                room_overlaps = []
                for (di, yr, sx, _, _) in all_sections:
                    for scx in course_subs[(di, yr, sx)]:
                        if not scx["is_lab_sub"]:
                            cd = scx["code"]
                            room_overlaps.append(theory_usage[(di, yr, sx, cd, d, t_i, r)])
                model.Add(sum(room_overlaps) <= 1)
            for lr in lab_rooms:
                lab_overlaps = []
                for (di, yr, sx, _, _) in all_sections:
                    for scx in course_subs[(di, yr, sx)]:
                        if scx["is_lab_sub"]:
                            cd = scx["code"]
                            lab_overlaps.append(lab_usage[(di, yr, sx, cd, d, t_i, lr)])
                model.Add(sum(lab_overlaps) <= 1)

    for (disc, yr, sec_idx, _, _) in all_sections:
        exs = course_subs[(disc, yr, sec_idx)]
        for d in days:
            for t_i in range(num_timeslots):
                model.Add(sum(schedule[(disc, yr, sec_idx, scx["code"], d, t_i)] for scx in exs) <= 1)

    def add_theory_constraints(model, sched, disc, yr, sec_idx, scode, needed):
        daySum = []
        for dd in days:
            dv = sum(sched[(disc, yr, sec_idx, scode, dd, t)] for t in range(num_timeslots))
            model.Add(dv <= 2)
            daySum.append(dv)
        model.Add(sum(daySum) == needed)

    def add_2consecutive(model, sched, disc, yr, sec_idx, scode):
        def block_of(t):
            if t < 3:
                return 1
            elif t < 5:
                return 2
            else:
                return 3
        allv = []
        for dd in days:
            for t_i in range(num_timeslots):
                allv.append(sched[(disc, yr, sec_idx, scode, dd, t_i)])
        model.Add(sum(allv) == 2)
        dayUsed = []
        for dd in days:
            ds = sum(sched[(disc, yr, sec_idx, scode, dd, t)] for t in range(num_timeslots))
            db = model.NewBoolVar(f"dayUsed_{disc}_{yr}_{sec_idx}_{scode}_{dd}")
            model.Add(ds == 2).OnlyEnforceIf(db)
            model.Add(ds == 0).OnlyEnforceIf(db.Not())
            dayUsed.append(db)
            pairVars = []
            for t_i in range(num_timeslots - 1):
                if block_of(t_i) == block_of(t_i+1):
                    p = model.NewBoolVar(f"pair_{disc}_{yr}_{sec_idx}_{scode}_{dd}_{t_i}")
                    model.Add(sched[(disc, yr, sec_idx, scode, dd, t_i)] == 1).OnlyEnforceIf(p)
                    model.Add(sched[(disc, yr, sec_idx, scode, dd, t_i+1)] == 1).OnlyEnforceIf(p)
                    for ot in range(num_timeslots):
                        if ot not in [t_i, t_i+1]:
                            model.Add(sched[(disc, yr, sec_idx, scode, dd, ot)] == 0).OnlyEnforceIf(p)
                    pairVars.append(p)
            model.Add(sum(pairVars) == 1).OnlyEnforceIf(db)
        model.Add(sum(dayUsed) == 1)

    for (disc, yr, sec_idx, _, _) in all_sections:
        exs = course_subs[(disc, yr, sec_idx)]
        base_map = defaultdict(list)
        for sc in exs:
            baseName = sc["code"].split("_")[0]
            base_map[baseName].append(sc["code"])
        for bName, codeList in base_map.items():
            t2_ = [x for x in codeList if "_T2Cons" in x]
            t1_ = [x for x in codeList if "_T1Sep" in x]
            if t2_ and t1_:
                t2c = t2_[0]
                t1c = t1_[0]
                for dd in days:
                    st2 = sum(schedule[(disc, yr, sec_idx, t2c, dd, ti)] for ti in range(num_timeslots))
                    st1 = sum(schedule[(disc, yr, sec_idx, t1c, dd, ti)] for ti in range(num_timeslots))
                    model.Add(st2 + st1 <= 2)

    for (disc, yr, sec_idx, _, _) in all_sections:
        exs = course_subs[(disc, yr, sec_idx)]
        for sc in exs:
            scode = sc["code"]
            is_lab = sc["is_lab_sub"]
            th = sc.get("theory_hours", 0)
            if is_lab and sc["needs_lab_2slots"]:
                add_2consecutive(model, schedule, disc, yr, sec_idx, scode)
            else:
                if th == 1:
                    model.Add(sum(schedule[(disc, yr, sec_idx, scode, dd, tt)]
                                  for dd in days for tt in range(num_timeslots)) == 1)
                elif th == 2 and sc["needs_lab_2slots"]:
                    add_2consecutive(model, schedule, disc, yr, sec_idx, scode)
                elif th in [3, 4]:
                    add_theory_constraints(model, schedule, disc, yr, sec_idx, scode, th)
                elif th == 2:
                    add_theory_constraints(model, schedule, disc, yr, sec_idx, scode, 2)

    teacher_prefs_info = []
    subcourse_map = {}
    for (disc, yr, sec_idx, _, _) in all_sections:
        for sc in course_subs[(disc, yr, sec_idx)]:
            subcourse_map[(disc, yr, sec_idx, sc["code"])] = (sc["code"].split("_")[0], sc["is_lab_sub"])

    for tID, pref_list in teacher_time_pref_list.items():
        for pref_dict in pref_list:
            # Now expect teacher prefs to use 'batch' as year and 'section' as section index.
            pref_year = pref_dict.get("batch", None)
            pref_section = pref_dict.get("section", None)
            c_name = pref_dict["course_name"]
            c_type = pref_dict["course_type"]
            day_pref = pref_dict.get("day", None)
            start_i = pref_dict.get("start_slot_index", None)
            end_i = pref_dict.get("end_slot_index", None)
            is_hard = pref_dict.get("is_hard", False)
            matched_subcourses = []
            for (disc, yr, sec_idx, scode, teacherCandidate) in teacher_assign:
                if teacherCandidate != tID:
                    continue
                if pref_year is not None:
                    if (disc, yr) != pref_year:
                        continue
                if pref_section is not None:
                    if sec_idx != pref_section:
                        continue
                baseN, isLab = subcourse_map[(disc, yr, sec_idx, scode)]
                if baseN == c_name and ((isLab and c_type == "lab") or (not isLab and c_type == "theory")):
                    matched_subcourses.append((disc, yr, sec_idx, scode))
            if is_hard and day_pref is not None and start_i is not None and end_i is not None:
                day_abbr = day_pref[:3]
                length_required = (end_i - start_i + 1)
                for (di, yr, sx, scd) in matched_subcourses:
                    assign_var = teacher_assign[(di, yr, sx, scd, tID)]
                    sc_len = subcourse_length(scd)
                    if sc_len == length_required:
                        for other_day in days:
                            if other_day != day_pref:
                                for t_i in range(num_timeslots):
                                    model.Add(schedule[(di, yr, sx, scd, other_day, t_i)] == 0).OnlyEnforceIf(assign_var)
                        for slot_i in range(num_timeslots):
                            if start_i <= slot_i <= end_i:
                                model.Add(schedule[(di, yr, sx, scd, day_pref, slot_i)] == 1).OnlyEnforceIf(assign_var)
                            else:
                                model.Add(schedule[(di, yr, sx, scd, day_pref, slot_i)] == 0).OnlyEnforceIf(assign_var)
            teacher_prefs_info.append((tID, pref_dict, matched_subcourses))

    for (disc, yr, sec_idx, _, _) in all_sections:
        for (disc, yr, sec_idx, scode, tID) in list(teacher_assign.keys()):
            ta_var = teacher_assign[(disc, yr, sec_idx, scode, tID)]
            pref_floor, floor_is_hard = teacher_info[tID]["floor_pref"]
            if floor_is_hard and pref_floor is not None:
                baseN, is_lab_sub = subcourse_map[(disc, yr, sec_idx, scode)]
                if is_lab_sub:
                    for lr in lab_rooms:
                        if lab_room_info[lr]["floor"] != pref_floor:
                            model.Add(lab_room_assign[(disc, yr, sec_idx, scode, lr)] == 0).OnlyEnforceIf(ta_var)
                else:
                    for d in days:
                        for t_i in range(num_timeslots):
                            for r in theory_rooms_original:
                                if room_info[r]["floor"] != pref_floor:
                                    model.Add(theory_usage[(disc, yr, sec_idx, scode, d, t_i, r)] == 0).OnlyEnforceIf(ta_var)

    hard_time_pref_subcourses = set()
    for tID, pref_dict, matched_subcourses in teacher_prefs_info:
        if pref_dict.get("is_hard", False) and \
          "day" in pref_dict and \
          "start_slot_index" in pref_dict and \
          "end_slot_index" in pref_dict:
            for (di, yr, sx, scd) in matched_subcourses:
                hard_time_pref_subcourses.add((di, yr, sx, scd))

    # Add day-of-week pattern constraints
    for (disc, year, section_index, batch_id, sec) in all_sections:
        if (disc, year, section_index) in skip_day_pref_for:
            continue

        exs = course_subs.get((disc, year, section_index), [])
        for sc in exs:
            scode = sc["code"]
            key = (disc, year, section_index, scode)

            # Skip if subcourse has hard time preference
            if key in hard_time_pref_subcourses:
                continue

            is_lab = sc["is_lab_sub"]
            year_val = year  # Assuming year is integer (1,2,3,4)

            if is_lab:
                if is_first_third_year(year_val):
                    # Labs for 1st/3rd years: Disallow Mon/Wed
                    for dd in ["Mon", "Wed"]:
                        for t_i in range(num_timeslots):
                            model.Add(
                                schedule[(disc, year, section_index, scode, dd, t_i)] == 0
                            )
                elif is_second_fourth_year(year_val):
                    # Labs for 2nd/4th years: Disallow Tue/Thu
                    for dd in ["Tue", "Thu"]:
                        for t_i in range(num_timeslots):
                            model.Add(
                                schedule[(disc, year, section_index, scode, dd, t_i)] == 0
                            )
            else:
                if is_first_third_year(year_val):
                    # Theory for 1st/3rd years: Disallow Tue/Thu
                    for dd in ["Tue", "Thu"]:
                        for t_i in range(num_timeslots):
                            model.Add(
                                schedule[(disc, year, section_index, scode, dd, t_i)] == 0
                            )
                elif is_second_fourth_year(year_val):
                    # Theory for 2nd/4th years: Disallow Mon/Wed
                    for dd in ["Mon", "Wed"]:
                        for t_i in range(num_timeslots):
                            model.Add(
                                schedule[(disc, year, section_index, scode, dd, t_i)] == 0
                            )
    for (disc, yr, sec_idx, _, _) in all_sections:
        if (disc, yr, sec_idx) not in skip_day_pref_for:
            model.Add(sum(theory_room_assign[(disc, yr, sec_idx, r)] for r in theory_rooms_original) == 1)

    for (disc, yr, sec_idx, _, _) in all_sections:
        for sc in course_subs[(disc, yr, sec_idx)]:
            if sc["is_lab_sub"]:
                scode = sc["code"]
                model.Add(sum(lab_room_assign[(disc, yr, sec_idx, scode, lr)] for lr in lab_rooms) == 1)
                for d in days:
                    for t_i in range(num_timeslots):
                        model.Add(sum(lab_usage[(disc, yr, sec_idx, scode, d, t_i, lr)] for lr in lab_rooms)
                                  == schedule[(disc, yr, sec_idx, scode, d, t_i)])
            else:
                scode = sc["code"]
                for d in days:
                    for t_i in range(num_timeslots):
                        model.Add(sum(theory_usage[(disc, yr, sec_idx, scode, d, t_i, r)] for r in theory_rooms_original)
                                  == schedule[(disc, yr, sec_idx, scode, d, t_i)])
    # --- Begin soft preference objective addition for soft (non-hard) preferences ---
    soft_penalties = []
    for (tID, pref_dict, matched_subcourses) in teacher_prefs_info:
        # Process only soft preferences with a defined day, start, and end.
        if (not pref_dict.get("is_hard", False) and
            pref_dict.get("day") is not None and
            pref_dict.get("start_slot_index") is not None and
            pref_dict.get("end_slot_index") is not None):
            day_pref = pref_dict.get("day")
            start_i = pref_dict.get("start_slot_index")
            end_i = pref_dict.get("end_slot_index")
            for (di, yr, sx, scd) in matched_subcourses:
                # Note: teacher_assign keys were stored as (disc, yr, sec_idx, scode, tid_str)
                assign_var = teacher_assign[(di, yr, sx, scd, str(tID))]
                for slot_i in range(start_i, end_i + 1):
                    sch_var = schedule[(di, yr, sx, scd, day_pref, slot_i)]
                    # Create an auxiliary integer variable (0 or 1) that is 1 if the teacher is assigned
                    # but the soft slot is NOT filled.
                    penalty = model.NewIntVar(0, 1, f"soft_penalty_{di}_{yr}_{sx}_{scd}_{tID}_{day_pref}_{slot_i}")
                    # The following constraints force penalty to equal 1 if (assign_var == 1 and sch_var == 0)
                    model.Add(penalty >= assign_var + (1 - sch_var) - 1)
                    model.Add(penalty <= assign_var)
                    model.Add(penalty <= 1 - sch_var)
                    soft_penalties.append(penalty)
    if soft_penalties:
        # Minimize the sum of penalties.
        model.Minimize(sum(soft_penalties))
    # --- End soft preference objective addition ---

    return (model, schedule, theory_room_assign, lab_room_assign, course_subs,
            all_sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info)

# 4) ADDITIONAL CONSTRAINT: MAX GAP
def add_max_gap_constraint(model, schedule, course_subs, all_sections, max_gap=1):
    for (disc, yr, sec_idx, _, _) in all_sections:
        exs = course_subs[(disc, yr, sec_idx)]
        for d in days:
            usedBools = []
            for t_i in range(num_timeslots):
                bvar = model.NewBoolVar(f"Used_{disc}_{yr}_{sec_idx}_{d}_{t_i}")
                ssum = sum(schedule[(disc, yr, sec_idx, sc["code"], d, t_i)] for sc in exs)
                model.Add(ssum >= bvar)
                model.Add(ssum <= len(exs) * bvar)
                usedBools.append(bvar)
            for t_i in range(num_timeslots - 1):
                anyLater = model.NewBoolVar(f"anyLater_{disc}_{yr}_{sec_idx}_{d}_{t_i}")
                model.Add(sum(usedBools[x] for x in range(t_i+1, num_timeslots)) >= 1).OnlyEnforceIf(anyLater)
                model.Add(sum(usedBools[x] for x in range(t_i+1, num_timeslots)) == 0).OnlyEnforceIf(anyLater.Not())
                end_lim = min(t_i + 1 + (max_gap + 1), num_timeslots)
                model.Add(sum(usedBools[x] for x in range(t_i+1, end_lim)) >= 1).OnlyEnforceIf([usedBools[t_i], anyLater])

# 5) BUILD AND SOLVE
def build_and_solve_once(disc_info, teacher_info, max_time_sec=60, skip_lastN_sections=0, max_gap=1, locked_slots=None, disabled_days=None):
    bscs1_count = disc_info.get(("BSCS", 1), {}).get("num_sections", 0)
    skip_day_pref_for = set()
    skip_single_room_for = set()
    if bscs1_count > 0 and skip_lastN_sections > 0:
        to_skip = []
        idx = bscs1_count - 1
        c = 0
        while c < skip_lastN_sections and idx >= 0:
            to_skip.append(idx)
            idx -= 1
            c += 1
        for i in to_skip:
            skip_day_pref_for.add(("BSCS", 1, i))
            skip_single_room_for.add(("BSCS", 1, i))
    (model, schedule, thrm, labrm, subs, sections,
     theory_usage, lab_usage, teacher_assign, teacher_prefs_info) = build_model_no_gap_with_teachers(
         disc_info, teacher_info, locked_slots=locked_slots, disabled_days=disabled_days,
         skip_day_pref_for=skip_day_pref_for,
         skip_single_room_for=skip_single_room_for
     )
    add_max_gap_constraint(model, schedule, subs, sections, max_gap)
    solver = cp_model.CpSolver()
    solver.parameters.log_search_progress = False
    solver.parameters.num_search_workers = 2
    solver.parameters.max_time_in_seconds = max_time_sec
    status = solver.Solve(model)
    feasible = (status in [cp_model.FEASIBLE, cp_model.OPTIMAL])
    return (feasible, solver, model, schedule, thrm, labrm, subs, sections,
            theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status)

def solve_iterative_with_timeouts(disc_info, teacher_info, locked_slots=None, disabled_days=None):
    if locked_slots is None:
        locked_slots = []
    if disabled_days is None:
        disabled_days = {}
    feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status = build_and_solve_once(
        disc_info, teacher_info, max_time_sec=60, skip_lastN_sections=0, max_gap=1,locked_slots=locked_slots, disabled_days=disabled_days )
    if feasible:
        return feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status

    feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status = build_and_solve_once(
        disc_info, teacher_info, max_time_sec=60, skip_lastN_sections=1, max_gap=1)
    if feasible:
        return feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status

    feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status = build_and_solve_once(
        disc_info, teacher_info, max_time_sec=60, skip_lastN_sections=2, max_gap=1)
    if feasible:
        return feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status

    feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status = build_and_solve_once(
        disc_info, teacher_info, max_time_sec=3000, skip_lastN_sections=3, max_gap=0)
    return feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status

# 6) EXTRACTION OF TIMETABLE INFORMATION
def extract_timetable(solver, schedule, thrm, labrm, subs, sections, teacher_assign, locked_slots_input, teacher_prefs_info, teacher_info):
    timetable_headers = []
    timetable_details = []
    detail_id = 4040
    timetable_id_counter = 3030
    day_full = {"Mon": "Monday", "Tue": "Tuesday", "Wed": "Wednesday", "Thu": "Thursday", "Fri": "Friday"}

    section_to_timetable_id = {}
    for (disc, year, sec_idx, batch_id, sec) in sections:
        t_id = timetable_id_counter
        timetable_id_counter += 1
        header = {
            "Timetable_ID": t_id,
            "Batch_ID": batch_id,  # for response
            "Section_ID": sec,     # for response
            "Status": "draft"
        }
        timetable_headers.append(header)
        section_to_timetable_id[(disc, year, sec_idx)] = t_id

    for (disc, year, sec_idx, batch_id, sec) in sections:
        t_id = section_to_timetable_id[(disc, year, sec_idx)]
        sub_list = subs[(disc, year, sec_idx)]
        for sub in sub_list:
            scode = sub["code"]
            teacher_id = None
            for key in teacher_assign:
                (d, yr, sx, sc, tID) = key
                if d == disc and yr == year and sx == sec_idx and sc == scode:
                    if solver.Value(teacher_assign[key]) > 0:
                        teacher_id = teacher_info[tID].get("id", tID)
                        break
            for day in days:
                active_slots = []
                for t_i in range(num_timeslots):
                    if solver.Value(schedule[(disc, year, sec_idx, scode, day, t_i)]) == 1:
                        active_slots.append(t_i)
                if not active_slots:
                    continue
                groups = []
                current_group = [active_slots[0]]
                for idx in active_slots[1:]:
                    if idx == current_group[-1] + 1:
                        current_group.append(idx)
                    else:
                        groups.append(current_group)
                        current_group = [idx]
                groups.append(current_group)
                for group in groups:
                    start_index = group[0]
                    end_index = group[-1]
                    start_time = timeslots[start_index].split('-')[0]
                    end_time = timeslots[end_index].split('-')[1]
                    if len(start_time) == 5:
                        start_time += ":00"
                    if len(end_time) == 5:
                        end_time += ":00"
                    room_id = None
                    if sub["is_lab_sub"]:
                        for lr in lab_rooms:
                            key_lr = (disc, year, sec_idx, scode, lr)
                            if key_lr in labrm and solver.Value(labrm[key_lr]) == 1:
                                room_id = lr
                                break
                    else:
                        for r in theory_rooms_original:
                            if (disc, year, sec_idx, r) in thrm and solver.Value(thrm[(disc, year, sec_idx, r)]) == 1:
                                room_id = r
                                break
                    locked = False
                    if locked_slots_input:
                        for lock in locked_slots_input:
                            if (lock.get("disc") == disc and lock.get("year") == year and
                                lock.get("section") == sec_idx and lock.get("scode") == scode and
                                lock.get("day") == day and lock.get("timeslot") in group):
                                locked = True
                                break
                    # pref_status = ""
                    if teacher_id is not None:
                      pref_status = ""  # Initialize as empty
                      for (tID, pref_dict, matched_subcourses) in teacher_prefs_info:
                          if tID != teacher_id or (disc, year, sec_idx, scode) not in matched_subcourses:
                              continue
                          if pref_dict.get("day") == day[:3]:
                              pref_start = pref_dict.get("start_slot_index")
                              pref_end = pref_dict.get("end_slot_index")
                              if pref_start is not None and pref_end is not None:
                                  pref_length = pref_end - pref_start + 1
                                  overlap = max(0, min(pref_end, end_index) - max(pref_start, start_index) + 1)
                                  if pref_dict.get("is_hard", False):
                                      pref_status = "H" if overlap == pref_length else ""
                                  else:
                                      if overlap == pref_length:
                                          pref_status = "G"
                                      elif overlap > 0:
                                          pref_status = "Y"
                                      else:
                                          pref_status = "R"
                              break  # Use the first matching teacher pref for the current day

                    detail = {
                        "Detail_ID": detail_id,
                        "Timetable_ID": t_id,
                        "Course_ID": sub.get("course_id", 1),
                        "Teacher_ID": teacher_id,
                        "Room_ID": room_id if room_id is not None else 1,
                        "Day": day_full.get(day, day),
                        "Start_Time": start_time,
                        "End_Time": end_time,
                        "Locked": locked,
                        "Teacher_pref_status": pref_status,
                        "Theory_or_Lab": "lab" if sub["is_lab_sub"] else "theory",
                        "Hard_slot": True if pref_status == "H" else False
                    }
                    detail_id += 1
                    timetable_details.append(detail)
    return timetable_headers, timetable_details

# 7) SETUP THE FLASK API SERVER AND NGROK TUNNEL
app = Flask(__name__)
CORS(app)
@app.route('/solve', methods=['POST'])
def solve_endpoint():
    data = request.get_json()
    room_info_input = data.get("rooms", {})
    lab_room_info_input = data.get("lab_rooms", {})
    teacher_info_input = data.get("teachers", {})
    discipline_info_input = data.get("discipline_info", {})
    locked_slots_input = data.get("locked_slots", [])
    # disabled_days_input = data.get("disabled_days", {})
    disabled_days_input = {}
    for key_str, days_list in data.get("disabled_days", {}).items():
        if key_str.startswith("(") and key_str.endswith(")"):
            # Corrected tuple key parsing
            clean_str = key_str[1:-1].strip()  # Remove parentheses and outer whitespace
            parts = [part.strip() for part in clean_str.split(',')]  # Keep inner spaces
            if len(parts) == 3:
                try:
                    discipline = parts[0]
                    year = int(parts[1])
                    section_idx = int(parts[2])
                    disabled_days_input[(discipline, year, section_idx)] = days_list
                except ValueError:
                    pass
    print(disabled_days_input)
    global room_info, lab_room_info, theory_rooms_original, lab_rooms
    room_info = room_info_input
    lab_room_info = lab_room_info_input
    theory_rooms_original = list(room_info.keys())
    lab_rooms = list(lab_room_info.keys())

    discipline_year_info = {}
    for key, value in discipline_info_input.items():
        if isinstance(key, str) and ',' in key:
            parts = key.split(',')
            disc = parts[0].strip()
            year = int(parts[1].strip())
            discipline_year_info[(disc, year)] = value
        else:
            discipline_year_info[key] = value

    feasible, solver, model, schedule, thrm, labrm, subs, sections, theory_usage, lab_usage, teacher_assign, teacher_prefs_info, status = solve_iterative_with_timeouts(discipline_year_info, teacher_info_input,locked_slots=locked_slots_input, disabled_days=disabled_days_input )

    if not feasible:
        return jsonify({"error": "No feasible solution found after all approaches."}), 400

    headers, details = extract_timetable(solver, schedule, thrm, labrm, subs, sections, teacher_assign, locked_slots_input, teacher_prefs_info, teacher_info_input)

    response = {
        "timetable_headers": headers,
        "timetable_details": details
    }
    return jsonify(response)

# if __name__ == '__main__':
#     ngrok.set_auth_token("2sxAj38no2f0yX7MeLxzZ7knBZh_4hf3Uehkva2sUMCX6UmWQ")
#     public_url = ngrok.connect(5000).public_url
#     print(" * ngrok tunnel:", public_url)
#     app.run(host="0.0.0.0", port=5000, debug=False, use_reloader=False)
if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000)

