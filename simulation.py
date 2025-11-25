import numpy as np
import heapq


# Hospital and Patient Parameters
lambda_1 = 10  # Arrival rate of patients (patients per time unit), exponential, mean is lambda_1 patients per time unit
lambda_2 = 5  # Arrival rate of patients (patients per time unit), exponential, mean is lambda_2 patients per time unit
mu1 = 1 / 0.598612  # Service rate in the ICU, lognormal, underlying normal mean is 1 / mu1
mu2 = 1 / 0.598612  # Service rate in the SDU, lognormal, underlying normal mean is 1 / mu1
sigma1 = 1  # Lognormal
sigma2 = 1  # Lognormal
c1 = 45  # Number of beds in the ICU
c2 = 35 # Number of beds in the SDU
x = 1.3  # Offplacement additional service time for critical
y = 1.2  # Offplacement additional service time for semi_critical
Q = 6  # Maximum critical queue length
K = Q + c1 + c2  # Maximum system capacity
theta = 1  # Abandonment rate

p = 0.30   # Probability of becoming semi_critical after leaving the critical state
p_cd = 0.05 # Probability of dying after leaving critical state
p_cd_offplace = 0.03 # Increase probability of dying after leaving critical state if offplaced
p_ccR = 0.07  # Probability of returning as critical after leaving the critical state
p_ccR_offplace = 0.05 # Increase probability of returning as critical after leaving critical state if offplaced
p_cscR = 0.07  # Probability of returning as semi_critical after leaving the critical state
p_cscR_offplace = 0.05 # Increase probability of returning as semi-critical after leaving critical state if offplaced

p_scd = 0.02 # Probability of dying after leaving semi-critical state
p_scd_offplace = 0.01 # Increased probability of dying after leaving semi-critical state if offplaced
p_scsc_R = 0.07  # Probability of returning as semi_critical after leaving the semi_critical state
p_scsc_R_offplace = 0.05 # Increase probability of returning as semi-critical after leaving semi-critical state if offplaced
p_scc_R = 0.07  # Probability of returning as critical after leaving the semi_critical state
p_scc_R_offplace = 0.05 # Increase probability of returning as critical after leaving semi-critical state if offplaced
delta = 1  # Return time, exponential, mean is delta

num_of_days = 5000  # Number of events to simulate (is not equal to the number of patients)


# Event queue (min-heap)
event_queue = []
heapq.heappush(event_queue, (0, "critical_arrival", 0))  # First critical arrival at time 0
heapq.heappush(event_queue,(0.1, "semi_critical_arrival", 1))  # First semi-critical arrival at time 0.1

# Tracking variables
critical_customers = set()
customers = {}
total_customers = 0
servers_1 = []
servers_2 = []
servers_3 = []
critical_queue = []
return_critical_queue = []
time_points = []
time_points_freq_1 = []
time_points_freq_2 = []
occupancy_1 = []
occupancy_2 = []
offplaced_occupancy_2 = []
offplaced_occupancy_3 = []
offplaced_occupancy_2_freq_1 = []
offplaced_occupancy_2_freq_2 = []
offplaced_occupancy_3_freq_1 = []
offplaced_occupancy_3_freq_2 = []
total_balked_times = 0
total_queue_times = 0
queue_times_count = 0
total_offplacement_times_critical = 0
total_offplacement_times_semi_critical = 0
offplacement_times_count_critical = 0
offplacement_times_count_semi_critical = 0
total_abandon_times = 0
total_bumping_times = 0
total_death_critical = 0
total_death_semi_critical = 0
total_return_critical = 0
total_return_semi_critical = 0
total_critical_treatments = 0
total_semi_critical_treatments = 0


# Start Simulation
event_time = 0
recording_flag = False
event_num = 0
record_time_1 = num_of_days / 10
record_time_2 = num_of_days / 10
while event_time <= num_of_days:
    print("")
    print(event_time)
    print(len(event_queue))
    print("")
    event_num += 1
    if event_num > 100000:
        break

    if event_time >= num_of_days / 10:
        recording_flag = True

    # Process next event
    event_time, event_type, customer_id = heapq.heappop(event_queue)

    while event_time >= record_time_1:
        offplaced_occupancy_2_freq_1.append(sum(1 for server_status in servers_2 if server_status.get("kind") == "1"))
        offplaced_occupancy_3_freq_1.append(len(servers_3))
        time_points_freq_1.append(record_time_1)
        record_time_1 += 1
    while event_time >= record_time_2:
        offplaced_occupancy_2_freq_2.append(sum(1 for server_status in servers_2 if server_status.get("kind") == "1"))
        offplaced_occupancy_3_freq_2.append(len(servers_3))
        time_points_freq_2.append(record_time_2)
        record_time_2 += 0.25
    time_points.append(event_time)
    occupancy_1.append(len(servers_1))
    occupancy_2.append(len(servers_2))
    offplaced_occupancy_2.append(sum(1 for server_status in servers_2 if server_status.get("kind") == "1"))
    offplaced_occupancy_3.append(len(servers_3))



    if event_type == "critical_arrival":
        customer_data = {
            "arrival_times": [event_time],
            "stage1_start_times": [],
            "stage1_service_times": [],
            "offplaced":[],
            "stage2_start_times": [],
            "stage2_service_times": [],
            "queue_times": [],
            "balked": False,
            "abandoned": False,
            "bumped_times": 0,
            "left": False,
            "died": False
        }
        print(f"new critical patient {customer_id} arrived.")

        if len(critical_queue) + len(return_critical_queue) >= Q:  # Queue length exceeded, balk
            customer_data["balked"] = True
            customers[customer_id] = customer_data
            if recording_flag:
                total_balked_times += 1
                total_customers += 1
                critical_customers.add(customer_id)
            next_arrival_time = event_time + np.random.exponential(1 / lambda_1)
            heapq.heappush(event_queue, (next_arrival_time, "critical_arrival", customer_id + 2))
            print(f"critical patient {customer_id} balked.")
            continue
        else:
            critical_queue.append(customer_id)
            customers[customer_id] = customer_data
            if recording_flag:
                total_customers += 1
                critical_customers.add(customer_id)
            abandonment_time = event_time + np.random.exponential(1 / theta)
            heapq.heappush(event_queue, (abandonment_time, "abandon", customer_id))

        if return_critical_queue:
            served_customer_id = return_critical_queue[0]
            return_critical_queue.remove(served_customer_id)
            return_patient = True
        else:
            served_customer_id = critical_queue[0]
            critical_queue.remove(served_customer_id)
            return_patient = False

        if len(servers_1) < c1:  # Free bed in ICU
            service_start_time = event_time
            offplaced = False
            print(f"critical patient {served_customer_id} begins treatment in ICU")
        elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) < c2):  # No free bed in ICU but there are people to bump to the SDU
            service_start_time = event_time
            offplaced = False
            print(f"critical patient {served_customer_id} begins treatment in ICU, bumping one to SDU.")
            for server_status in servers_1[:]:
                if server_status["kind"] == "2":
                    bumped_customer_id = server_status["customer_id"]
                    heapq.heappush(event_queue,(event_time, "bump_to_SDU", bumped_customer_id))
                    servers_1.remove(server_status)
                    break
        elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) >= c2):  # No free bed in ICU or SDU but there are people to bump to the ward
            service_start_time = event_time
            offplaced = False
            print(f"critical patient {served_customer_id} begins treatment in ICU, bumping one to ward.")
            remaining_time = float('inf')
            for server_status in servers_1[:]:
                if (server_status["kind"] == "2") & (server_status["time"] < remaining_time):
                    bumped_customer_id = server_status["customer_id"]
                    server_status_to_be_removed = server_status
            heapq.heappush(event_queue,(event_time, "bump_to_ward", bumped_customer_id))
            servers_1.remove(server_status_to_be_removed)
            if recording_flag:
                total_bumping_times += 1
        elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) >= c1) and (len(servers_2) < c2):  # No free bed in ICU but free bed in SDU, offplace
            service_start_time = event_time
            offplaced = True
            print(f"critical patient {served_customer_id} begins treatment in SDU, offplaced.")
        else:  # No way to begin treatment, put back in queue
            if return_patient:
                return_critical_queue.insert(0, served_customer_id)
                next_arrival_time = event_time + np.random.exponential(1 / lambda_1)
                heapq.heappush(event_queue, (next_arrival_time, "critical_arrival", customer_id + 2))
                continue
            else:
                critical_queue.insert(0, served_customer_id)
                next_arrival_time = event_time + np.random.exponential(1 / lambda_1)
                heapq.heappush(event_queue, (next_arrival_time, "critical_arrival", customer_id + 2))
                continue

        if not offplaced :  # Patient is not balked nor offplaced
            queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
            service_time1 = np.random.lognormal(1/mu1, sigma1)
            customers[served_customer_id]["stage1_start_times"].append(service_start_time)
            customers[served_customer_id]["queue_times"].append(queue_time)
            customers[served_customer_id]["stage1_service_times"].append(service_time1)
            customers[served_customer_id]["offplaced"].append(offplaced)
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", served_customer_id))
            servers_1.append({"time": leave_time, "kind": "1", "customer_id": served_customer_id})
            print(f"updated patient {served_customer_id} information")

        else:  # Patient is offplaced
            print("OFFPLACE TIME ADDED_1")
            queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
            service_time1 = np.random.lognormal(1/mu1, sigma1) * x
            customers[served_customer_id]["stage1_start_times"].append(service_start_time)
            customers[served_customer_id]["queue_times"].append(queue_time)
            customers[served_customer_id]["stage1_service_times"].append(service_time1)
            customers[served_customer_id]["offplaced"].append(offplaced)
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
                total_offplacement_times_critical += service_time1
                offplacement_times_count_critical += 1
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", served_customer_id))
            servers_2.append({"time": leave_time, "kind": "1", "customer_id": served_customer_id})
            print(f"updated patient {served_customer_id} information")


        # Schedule next arrival
        next_arrival_time = event_time + np.random.exponential(1 / lambda_1)
        heapq.heappush(event_queue, (next_arrival_time, "critical_arrival", customer_id + 2))




    elif event_type == "return_critical":
        if recording_flag:
            total_return_critical += 1
        customers[customer_id]["arrival_times"].append(event_time)
        print(f"returning critical patient {customer_id} arrived.")

        return_critical_queue.append(customer_id)
        abandonment_time = event_time + np.random.exponential(1 / theta)
        critical_customers.add(customer_id)
        heapq.heappush(event_queue, (abandonment_time, "abandon", customer_id))
        served_customer_id = return_critical_queue[0]
        return_critical_queue.remove(served_customer_id)


        if len(servers_1) < c1:  # Free bed in ICU
            service_start_time = event_time
            offplaced = False
            print(f"critical patient {served_customer_id} begins treatment in ICU")
        elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) < c2):  # No free bed in ICU but there are people to bump to the SDU
            service_start_time = event_time
            offplaced = False
            print(f"critical patient {served_customer_id} begins treatment in ICU, bumping one to SDU")
            for server_status in servers_1[:]:
                if server_status["kind"] == "2":
                    bumped_customer_id = server_status["customer_id"]
                    heapq.heappush(event_queue,(event_time, "bump_to_SDU", bumped_customer_id))
                    servers_1.remove(server_status)
                    break
        elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) >= c2):  # No free bed in ICU or SDU but there are people to bump to the ward
            service_start_time = event_time
            offplaced = False
            print(f"critical patient {served_customer_id} begins treatment in ICU, bumping one to ward")
            remaining_time = float('inf')
            for server_status in servers_1[:]:
                if (server_status["kind"] == "2") & (server_status["time"] < remaining_time):
                    bumped_customer_id = server_status["customer_id"]
                    server_status_to_be_removed = server_status
            heapq.heappush(event_queue,(event_time, "bump_to_ward", bumped_customer_id))
            servers_1.remove(server_status_to_be_removed)
            if recording_flag:
                total_bumping_times += 1
        elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) == c1) and (len(servers_2) < c2):  # No free bed in ICU but free bed in SDU, offplace
            service_start_time = event_time
            offplaced = True
            print(f"critical patient begins treatment in SDU, offplaced")
        else:  # No way to begin treatment, put back to queue
            return_critical_queue.insert(0, served_customer_id)
            continue

        if not offplaced :  # Patient is not balked nor offplaced
            queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
            service_time1 = np.random.lognormal(1/mu1, sigma1)
            customers[served_customer_id]["stage1_start_times"].append(service_start_time)
            customers[served_customer_id]["queue_times"].append(queue_time)
            customers[served_customer_id]["stage1_service_times"].append(service_time1)
            customers[served_customer_id]["offplaced"].append(offplaced)
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", served_customer_id))
            servers_1.append({"time": leave_time, "kind": "1", "customer_id": served_customer_id})
            print(f"updated patient {served_customer_id} information")

        else:  # Patient is offplaced
            print("OFFPLACEMENT TIME ADDED_2")
            queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
            service_time1 = np.random.lognormal(1/mu1, sigma1) * y
            customers[served_customer_id]["stage1_start_times"].append(service_start_time)
            customers[served_customer_id]["queue_times"].append(queue_time)
            customers[served_customer_id]["stage1_service_times"].append(service_time1)
            customers[served_customer_id]["offplaced"].append(offplaced)
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
                total_offplacement_times_critical += service_time1
                offplacement_times_count_critical += 1
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", served_customer_id))
            servers_2.append({"time": leave_time, "kind": "1", "customer_id": served_customer_id})
            print(f"updated patient {served_customer_id} information")




    elif event_type == "semi_critical_arrival":
        customer_data = {
            "arrival_times": [event_time],
            "stage1_start_times": [],
            "stage1_service_times": [],
            "offplaced":[],
            "stage2_start_times": [],
            "stage2_service_times": [],
            "queue_times": [],
            "balked": False,
            "abandoned": False,
            "bumped_times": 0,
            "left": False,
            "died": False
        }

        customers[customer_id] = customer_data
        served_customer_id = customer_id
        if recording_flag:
            total_customers += 1
        print(f"new semi_critical patient {customer_id} arrived.")

        if len(servers_2) < c2:  # Free bed in SDU
            print(servers_2)
            print(len(servers_2) < c2)
            print(c2)
            print(len(servers_2))
            service_start_time = event_time
            offplaced = False
            service_time2 = np.random.lognormal(1/mu2, sigma2)
            customers[served_customer_id]["stage2_start_times"].append(service_start_time)
            customers[served_customer_id]["stage2_service_times"].append(service_time2)
            customers[served_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time2
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
            servers_2.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            print(f"semi_critical patient {served_customer_id} begins treatment in SDU")
            print(f"updated patient {served_customer_id} information")
        elif len(servers_1) < c1:  # No Free bed in SDU, free bed in ICU
            service_start_time = event_time
            offplaced = False
            service_time2 = np.random.lognormal(1/mu2, sigma2)
            customers[served_customer_id]["stage2_start_times"].append(service_start_time)
            customers[served_customer_id]["stage2_service_times"].append(service_time2)
            customers[served_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time2
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
            servers_1.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            print(f"semi_critical patient {served_customer_id} begins treatment in ICU")
            print(f"updated patient {served_customer_id} information")
        elif len(servers_1) >= c1 and len(servers_2) >= c2:  # No free bed in ICU and SDU, offplace
            service_start_time = event_time
            offplaced = True
            service_time2 = np.random.lognormal(1/mu2, sigma2)
            customers[served_customer_id]["stage2_start_times"].append(service_start_time)
            customers[served_customer_id]["stage2_service_times"].append(service_time2)
            customers[served_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time2
            if recording_flag:
                total_offplacement_times_semi_critical += service_time2
                offplacement_times_count_semi_critical += 1
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
            servers_3.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            print(f"semi_critical patient {served_customer_id} begins treatment in ward, offplaced")
            print(f"updated patient {served_customer_id} information")

        # Schedule next arrival
        next_arrival_time = event_time + np.random.exponential(1 / lambda_2)
        heapq.heappush(event_queue, (next_arrival_time, "semi_critical_arrival", customer_id + 2))




    elif event_type == "return_semi_critical":
        if recording_flag:
            total_return_semi_critical += 1
        customers[customer_id]["arrival_times"].append(event_time)
        served_customer_id = customer_id
        print(f"returning semi_critical patient {customer_id} arrived.")

        if len(servers_2) < c2:  # Free bed in SDU
            service_start_time = event_time
            offplaced = False
            service_time2 = np.random.lognormal(1/mu2, sigma2)
            customers[served_customer_id]["stage2_start_times"].append(service_start_time)
            customers[served_customer_id]["stage2_service_times"].append(service_time2)
            customers[served_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time2
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
            servers_2.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            print(f"semi_critical patient {served_customer_id} begins treatment in SDU")
            print(f"updated patient {served_customer_id} information")
        elif len(servers_1) < c1:  # No free bed in SDU, free bed in ICU
            service_start_time = event_time
            offplaced = False
            service_time2 = np.random.lognormal(1/mu2, sigma2)
            customers[served_customer_id]["stage2_start_times"].append(service_start_time)
            customers[served_customer_id]["stage2_service_times"].append(service_time2)
            customers[served_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time2
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
            servers_1.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            print(f"semi_critical patient {served_customer_id} begins treatment in ICU")
            print(f"updated patient {served_customer_id} information")
        elif len(servers_1) >= c1 and len(servers_2) >= c2:  # No free bed in ICU and SDU, offplace
            service_start_time = event_time
            offplaced = True
            service_time2 = np.random.lognormal(1/mu2, sigma2) * y
            customers[served_customer_id]["stage2_start_times"].append(service_start_time)
            customers[served_customer_id]["stage2_service_times"].append(service_time2)
            customers[served_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time2
            if recording_flag:
                total_offplacement_times_semi_critical += service_time2
                offplacement_times_count_semi_critical += 1
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
            servers_3.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            print(f"semi_critical patient {served_customer_id} begins treatment in ward, offplaced")
            print(f"updated patient {served_customer_id} information")




    elif event_type == "critical_to_semi_critical":
        print(f"critical patient {customer_id} transitions to semi_critical")
        service_time2 = np.random.lognormal(1/mu2, sigma2)
        leave_time = service_time2 + event_time
        offplaced = False
        customers[customer_id]["stage2_start_times"].append(event_time)
        customers[customer_id]["stage2_service_times"].append(service_time2)
        customers[customer_id]["offplaced"].append(offplaced)
        for server_status in servers_1[:]:
                if (server_status["customer_id"] == customer_id) and (server_status["kind"] == "1"):
                    servers_1.remove(server_status)
                    servers_1.append({"time": leave_time, "kind": "2", "customer_id": customer_id})
                    print(f"{customer_id} is in ICU")
                    break
        for server_status in servers_2[:]:
                if (server_status["customer_id"] == customer_id) and (server_status["kind"] == "1"):
                    servers_2.remove(server_status)
                    servers_2.append({"time": leave_time, "kind": "2", "customer_id": customer_id})
                    print(f"{customer_id} is in SDU")
                    break
        heapq.heappush(event_queue, (leave_time, "end_semi_critical", customer_id))




    elif event_type == "abandon":
        if customer_id in critical_queue:
            critical_queue.remove(customer_id)
            customers[customer_id]["abandoned"] = True
            queue_time = event_time - customers[customer_id]["arrival_times"][-1]
            customers[customer_id]["queue_times"].append(queue_time)
            if recording_flag:
                total_abandon_times += 1
                total_queue_times += queue_time
                queue_times_count += 1
        if customer_id in return_critical_queue:
            return_critical_queue.remove(customer_id)
            customers[customer_id]["abandoned"] = True
            queue_time = event_time - customers[customer_id]["arrival_times"][-1]
            customers[customer_id]["queue_times"].append(queue_time)
            if recording_flag:
                total_abandon_times += 1
                total_queue_times += queue_time
                queue_times_count += 1




    elif event_type == "bump_to_SDU":
        print(f"semi_critical patient {customer_id} is bumped to the SDU")
        leave_time = customers[customer_id]["stage2_start_times"][-1] + customers[customer_id]["stage2_service_times"][-1]
        servers_2.append({"time": leave_time, "kind": 2, "customer_id": customer_id})




    elif event_type == "bump_to_ward":
        print(f"semi_critical patient {customer_id} is bumped to the ward")
        remaining_time = customers[customer_id]["stage2_start_times"][-1] + customers[customer_id]["stage2_service_times"][-1] - event_time
        leave_time = event_time + remaining_time * y
        servers_3.append({"time": leave_time, "kind": 2, "customer_id": customer_id})
        customers[customer_id]["bumped_times"] += 1
        customers[customer_id]["stage2_service_times"][-1] += remaining_time * (y - 1)
        if recording_flag:
            total_offplacement_times_semi_critical += remaining_time * y
            offplacement_times_count_semi_critical += 1
        heapq.heappush(event_queue, (event_time + remaining_time * y, "end_semi_critical", customer_id))




    elif event_type == "leave_ICU":
        customer_found = False
        for server_status in servers_1[:]:
            if server_status["customer_id"] == customer_id:
                servers_1.remove(server_status)
                customer_found = True
                print(f"removed patient {customer_id} from ICU")
                break
        if not customer_found:
            continue
        if sum(1 for server_status in servers_2 if (server_status["kind"] == "1")) > 0:  # Transfer offplaced critical patient to ICU
            longest_remaining_time = 0
            for server_status in servers_2:
                if server_status["kind"] == "1":
                    remaining_time = customers[server_status["customer_id"]]["stage1_start_times"][-1] + customers[server_status["customer_id"]]["stage1_service_times"][-1] - event_time
                    if remaining_time > longest_remaining_time:
                        longest_remaining_time = remaining_time
                        transfer_customer_id = server_status["customer_id"]
            for server_status in servers_2[:]:
                if server_status["customer_id"] == transfer_customer_id:
                    servers_2.remove(server_status)
                    break
            servers_1.append({"time": event_time + longest_remaining_time / x, "kind": "1", "customer_id": transfer_customer_id})
            print(f"critical patient {transfer_customer_id} is tranferred into the ICU")
            customers[transfer_customer_id]["stage1_service_times"][-1] = event_time + longest_remaining_time / x - customers[transfer_customer_id]["stage1_start_times"][-1]
            print("OFFPLACEMENT TIME ADDED_3")
            if recording_flag:
                total_offplacement_times_critical -= longest_remaining_time
            heapq.heappush(event_queue, (event_time + longest_remaining_time / x, "end_critical", transfer_customer_id))
        elif return_critical_queue:  # Returned critical patient in the queue, priority
            new_customer_id = return_critical_queue[0]
            return_critical_queue.remove(new_customer_id)
            service_start_time = event_time
            queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
            service_time1 = np.random.lognormal(1/mu1, sigma1)
            offplaced = False
            customers[new_customer_id]["stage1_start_times"].append(service_start_time)
            customers[new_customer_id]["queue_times"].append(queue_time)
            customers[new_customer_id]["stage1_service_times"].append(service_time1)
            customers[new_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", new_customer_id))
            servers_1.append({"time": leave_time, "kind": "1", "customer_id": new_customer_id})
            print(f"critical patient {new_customer_id} begins treatment in ICU")
        elif critical_queue:  # New critical patient in the queue
            new_customer_id = critical_queue[0]
            critical_queue.remove(new_customer_id)
            service_start_time = event_time
            queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
            service_time1 = np.random.lognormal(1/mu1, sigma1)
            offplaced = False
            customers[new_customer_id]["stage1_start_times"].append(service_start_time)
            customers[new_customer_id]["queue_times"].append(queue_time)
            customers[new_customer_id]["stage1_service_times"].append(service_time1)
            customers[new_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", new_customer_id))
            servers_1.append({"time": leave_time, "kind": "1", "customer_id": new_customer_id})
            print(f"critical patient {served_customer_id} begins treatment in ICU")
        elif servers_3:  # Transfer offplaced semi-critical patient to the ICU
            longest_remaining_time = 0
            for server_status in servers_3[:]:
                remaining_time = customers[server_status["customer_id"]]["stage2_start_times"][-1] + customers[server_status["customer_id"]]["stage2_service_times"][-1] - event_time
                if remaining_time > longest_remaining_time:
                    longest_remaining_time = remaining_time
                    transfer_customer_id = server_status["customer_id"]
            for server_status in servers_3[:]:
                if server_status["customer_id"] == transfer_customer_id:
                    servers_3.remove(server_status)
                    break
            if recording_flag:
                total_offplacement_times_semi_critical -= longest_remaining_time
            servers_1.append({"time": event_time + longest_remaining_time / y, "kind": "2", "customer_id": transfer_customer_id})
            customers[transfer_customer_id]["stage2_service_times"][-1] = event_time + longest_remaining_time / y - customers[transfer_customer_id]["stage2_start_times"][-1]
            heapq.heappush(event_queue, (event_time + longest_remaining_time / y, "end_semi_critical", transfer_customer_id))
            print(f"semi_critical patient {transfer_customer_id} is transferred into the ICU")




    elif event_type == "leave_SDU":
        customer_found = False
        for server_status in servers_2[:]:
            if server_status["customer_id"] == customer_id:
                servers_2.remove(server_status)
                print(f"removed patient {customer_id} from SDU")
                customer_found = True
                break
        if not customer_found:
            continue
        if servers_3:  # Transfer offplaced semi_critical patient to SDU
            longest_remaining_time = 0
            print(f"    ward: {servers_3}")
            for server_status in servers_3[:]:
                remaining_time = customers[server_status["customer_id"]]["stage2_start_times"][-1] + customers[server_status["customer_id"]]["stage2_service_times"][-1] - event_time
                if remaining_time > longest_remaining_time:
                    longest_remaining_time = remaining_time
                    transfer_customer_id = server_status["customer_id"]
            for server_status in servers_3[:]:
                if server_status["customer_id"] == transfer_customer_id:
                    servers_3.remove(server_status)
                    break
            if recording_flag:
                total_offplacement_times_semi_critical -= longest_remaining_time
            servers_2.append({"time": event_time + longest_remaining_time / y, "kind": "2", "customer_id": transfer_customer_id})
            customers[transfer_customer_id]["stage2_service_times"][-1] = event_time + longest_remaining_time / y - customers[transfer_customer_id]["stage2_start_times"][-1]
            heapq.heappush(event_queue, (event_time + longest_remaining_time / y, "end_semi_critical", transfer_customer_id))
            print(f"semi_critical patient {transfer_customer_id} is tranferred into the SDU")
            print(f"    ward: {servers_3}")
        elif return_critical_queue:  # Offplace returned critical patients
            new_customer_id = return_critical_queue[0]
            return_critical_queue.remove(new_customer_id)
            service_start_time = event_time
            queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
            service_time1 = np.random.lognormal(1/mu1, sigma1) * x
            print("OFFPLACEMENT TIME ADDED_4")
            if recording_flag:
                total_offplacement_times_critical += service_time1
                offplacement_times_count_critical += 1
            offplaced = True
            customers[new_customer_id]["stage1_start_times"].append(service_start_time)
            customers[new_customer_id]["queue_times"].append(queue_time)
            customers[new_customer_id]["stage1_service_times"].append(service_time1)
            customers[new_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", new_customer_id))
            servers_2.append({"time": leave_time, "kind": "1", "customer_id": new_customer_id})
            print(f"critical patient {new_customer_id} begins treatment in SDU, offplaced.")
            print("8")
        elif critical_queue:  # Offplace new critical patients
            new_customer_id = critical_queue[0]
            critical_queue.remove(new_customer_id)
            service_start_time = event_time
            queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
            if recording_flag:
                total_queue_times += queue_time
                queue_times_count += 1
            service_time1 = np.random.lognormal(1/mu1, sigma1) * x
            print("OFFPLACEMENT TIME ADDED_5")
            if recording_flag:
                total_offplacement_times_critical += service_time1
                offplacement_times_count_critical += 1
            offplaced = True
            customers[new_customer_id]["stage1_start_times"].append(service_start_time)
            customers[new_customer_id]["queue_times"].append(queue_time)
            customers[new_customer_id]["stage1_service_times"].append(service_time1)
            customers[new_customer_id]["offplaced"].append(offplaced)
            leave_time = service_start_time + service_time1
            heapq.heappush(event_queue, (leave_time, "end_critical", new_customer_id))
            servers_2.append({"time": leave_time, "kind": "1", "customer_id": new_customer_id})
            print(f"critical patient {new_customer_id} begins treatment in SDU, offplaced.")
            print("9")




    elif event_type == "leave_ward":
        print(f"Attempting to remove patient {customer_id} from ward")
        customer_found = False
        for server_status in servers_3[:]:
            if server_status["customer_id"] == customer_id:
                servers_3.remove(server_status)
                customer_found = True
                print(f"removed patient {customer_id} from ward")
                if recording_flag:
                    total_bumping_times += 1
                break




    elif event_type == "end_critical":
        print(f"checking if patient {customer_id} ended critical treatment")
        end_time_1 = 0
        end_time_2 = 0
        try:
            end_time_1 = customers[customer_id]["stage1_start_times"][-1] + customers[customer_id]["stage1_service_times"][-1]

        except Exception as e:
            pass

        if (abs(end_time_1 - event_time) < 1e-10) or (abs(end_time_2 - event_time) < 1e-10):
            decision = np.random.rand()
            print(f"patient {customer_id} ended critical treatment")
            offplace_penalty = customers[customer_id]["offplaced"][-1]
            if decision < p:
                # Schedule transition to semi-critical
                add_returning_event_flag = True
                for event in event_queue:
                    if (event[2] == customer_id) and ((event[1] == "return_critical") or (event[1] == "return_semi_critical") or (event[1] == "critical_to_semi_critical")):
                        add_returning_event_flag = False
                if add_returning_event_flag:
                    heapq.heappush(event_queue, (event_time, "critical_to_semi_critical", customer_id))
            elif decision < p + p_ccR + p_ccR_offplace * offplace_penalty:
                # Schedule leave and return as critical
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                add_returning_event_flag = True
                for event in event_queue:
                    if (event[2] == customer_id) and ((event[1] == "return_critical") or (event[1] == "return_semi_critical") or (event[1] == "critical_to_semi_critical")):
                        add_returning_event_flag = False
                if add_returning_event_flag:
                    return_time = event_time + np.random.exponential(delta)
                    heapq.heappush(event_queue, (return_time, "return_critical", customer_id))
            elif decision < p + p_ccR + p_ccR_offplace * offplace_penalty + p_cscR + p_cscR_offplace * offplace_penalty:
                # Schedule leave and return as semi-critical
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                add_returning_event_flag = True
                for event in event_queue:
                    if (event[2] == customer_id) and ((event[1] == "return_critical") or (event[1] == "return_semi_critical") or (event[1] == "critical_to_semi_critical")):
                        add_returning_event_flag = False
                if add_returning_event_flag:
                    return_time = event_time + np.random.exponential(delta)
                    heapq.heappush(event_queue, (return_time, "return_semi_critical", customer_id))
            elif decision < p + p_ccR + p_ccR_offplace * offplace_penalty + p_cscR + p_cscR_offplace * offplace_penalty + p_cd + p_cd_offplace * offplace_penalty:
                customers[customer_id]["died"] = True
                if recording_flag:
                    total_death_critical += 1
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
            else:
                # Schedule leave without returning to the system
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))




    elif event_type == "end_semi_critical":
        print(f"checking if patient {customer_id} ended semi_critical treatment")
        end_time_1 = 0
        end_time_2 = 0
        try:
            end_time_2 = customers[customer_id]["stage2_start_times"][-1] + customers[customer_id]["stage2_service_times"][-1]

        except Exception as e:
            pass

        if (abs(end_time_1 - event_time) < 1e-10) or (abs(end_time_2 - event_time) < 1e-10):
            decision = np.random.rand()
            print(f"patient {customer_id} ended semi_critical treatment")
            offplace_penalty = customers[customer_id]["offplaced"][-1]
            if decision < p_scsc_R + p_scsc_R_offplace * offplace_penalty:
                # Schedule leave and return as semi-critical
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ward", customer_id))
                add_returning_event_flag = True
                for event in event_queue:
                    if (event[2] == customer_id) and ((event[1] == "return_critical") or (event[1] == "return_semi_critical") or (event[1] == "critical_to_semi_critical")):
                        add_returning_event_flag = False
                if add_returning_event_flag:
                    return_time = event_time + np.random.exponential(delta)
                    heapq.heappush(event_queue, (event_time, "return_semi_critical", customer_id))
            elif decision < p_scsc_R + p_scsc_R_offplace * offplace_penalty + p_scc_R + p_scc_R_offplace * offplace_penalty:
                # Schedule leave and return as critical
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ward", customer_id))
                add_returning_event_flag = True
                for event in event_queue:
                    if (event[2] == customer_id) and ((event[1] == "return_critical") or (event[1] == "return_semi_critical") or (event[1] == "critical_to_semi_critical")):
                        add_returning_event_flag = False
                if add_returning_event_flag:
                    return_time = event_time + np.random.exponential(delta)
                    heapq.heappush(event_queue, (event_time, "return_critical", customer_id))
            elif decision < p_scsc_R + p_scsc_R_offplace * offplace_penalty + p_scc_R + p_scc_R_offplace * offplace_penalty + p_scd + p_scd_offplace * offplace_penalty:
                customers[customer_id]["died"] = True
                if recording_flag:
                    total_death_semi_critical += 1
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ward", customer_id))
            else:
                # Schedule leave without returning to the system
                heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                heapq.heappush(event_queue,(event_time, "leave_ward", customer_id))




# Statistics
total_customers = len(customers)

total_critical_customers = len(critical_customers)

balked_customers = sum(1 for c in customers.values() if c["balked"])

abandoned_customers = sum(1 for c in customers.values() if c["abandoned"])

bumped_customers = sum(1 for c in customers.values() if c["bumped_times"])

total_queue_time = []
for c in customers:
    for queue_time in customers[c]["queue_times"]:
        total_queue_time.append(queue_time)

count = 0
total_stage1_time = 0
for c in customers:
    for stage1_time in customers[c]["stage1_service_times"]:
        count += 1
        total_stage1_time += stage1_time
average_stage1_time = total_stage1_time / count

count = 0
total_stage2_time = 0
for c in customers:
    for stage2_time in customers[c]["stage2_service_times"]:
        count += 1
        total_stage2_time += stage2_time
average_stage2_time = total_stage2_time / count


print("")
print("----------")
print(f"Total patients: {total_customers}")
print(f"Total balked patients: {total_balked_times} ({balked_customers / total_critical_customers:.2%})")
print(f"Total abandoned patients: {abandoned_customers} ({abandoned_customers / total_critical_customers:.2%})")
print(f"Average queue time: {np.mean(total_queue_time):.4f}")
print(f"Average service time (ICU): {average_stage1_time:.4f}")
print(f"Average service time (SDU): {average_stage2_time:.4f}")
print(f"Average Number of Offplaced Critical Patients (1 day sampling): {np.mean(offplaced_occupancy_2_freq_1):.4f}")
print(f"Average Number of Offplaced Critical Patients (6 hour sampling): {np.mean(offplaced_occupancy_2_freq_2):.4f}")
print(f"Average Number of Offplaced Semi-Critical Patients (1 day sampling): {np.mean(offplaced_occupancy_3_freq_1):.4f}")
print(f"Average Number of Offplaced Semi-Critical Patients (6 hour sampling): {np.mean(offplaced_occupancy_3_freq_2):.4f}")
print(f"Total returning critical patients: {total_return_critical} ({total_return_critical / total_customers:.2%})")
print(f"Total returning semi-critical patients: {total_return_semi_critical} ({total_return_semi_critical / total_customers:.2%})")
print(f"Total died critical patients: {total_death_critical} ({total_death_critical / total_customers:.2%})")
print(f"Total died semi-critical patients: {total_death_semi_critical} ({total_death_semi_critical / total_customers:.2%})")


print(len(servers_2))
print(len(servers_3))
print("---------")
print("")


# Plots
import matplotlib.pyplot as plt

plt.hist([queue_time for queue_time in total_queue_time if queue_time >= 1e-3], bins=30, alpha=0.7, edgecolor='black')
plt.xlabel("Queue Time")
plt.ylabel("Frequency")
plt.title("Histogram of Queue Times")
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points, occupancy_1, label="Number of Patients in ICU", color='red', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("ICU Occupancy")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points, occupancy_2, label="Number of Patients in SDU", color='blue', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("SDU Occupancy")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points, offplaced_occupancy_2, label="Number of Offplaced Critical Patients", color='red', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("Offplaced Critical Patients")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points_freq_1, offplaced_occupancy_2_freq_1, label="Number of Offplaced Critical Patients", color='red', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("Offplaced Critical Patients (1 day sampling)")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points_freq_2, offplaced_occupancy_2_freq_2, label="Number of Offplaced Critical Patients", color='red', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("Offplaced Critical Patients (6 hour sampling)")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points, offplaced_occupancy_3, label="Number of Offplaced Semi-Critical Patients", color='blue', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("Offplaced Semi-Critical Patients")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points_freq_1, offplaced_occupancy_3_freq_1, label="Number of Offplaced Semi-Critical Patients", color='blue', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("Offplaced Semi-Critical Patients (1 day sampling)")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 5))
plt.scatter(time_points_freq_2, offplaced_occupancy_3_freq_2, label="Number of Offplaced Semi-Critical Patients", color='blue', alpha=0.6, s=3)
plt.xlabel("Day")
plt.ylabel("Number of Patients")
plt.title("Offplaced Semi-Critical Patients (6 hour sampling)")
plt.legend()
plt.grid()
plt.show()



# Store all patients' information for further inspection  
import pickle
with open("patient_information.pkl", "wb") as f:
    pickle.dump(customers, f)
