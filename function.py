import numpy as np
import heapq

def run_simulation(lambda_1 = 10.0, 
                   lambda_2 = 5.0, 
                   mu1 = 0.598612, 
                   mu2 = 0.598612, 
                   sigma1 = 1.0, 
                   sigma2 = 1.0, 
                   c1 = 23, 
                   c2 = 30, 
                   Q = 5, 
                   theta = 1.0,
                   p = 0.30, 
                   p_cd = 0.20, 
                   p_cd_offplace = 0.10, 
                   p_cc_R = 0.07, 
                   p_cc_R_offplace = 0.05, 
                   p_csc_R = 0.07, 
                   p_csc_R_offplace = 0.05, 
                   p_scd = 0.10, 
                   p_scd_offplace = 0.05, 
                   p_scc_R = 0.07, 
                   p_scc_R_offplace = 0.05, 
                   p_scsc_R = 0.07, 
                   p_scsc_R_offplace = 0.05, 
                   delta = 1.0, 
                   x = 1.3, 
                   y = 1.2, 
                   num_of_days = 10000, 
                   sampling_freq = 1.0):
    """
Run the simulation and return data about the simulation.


Parameters\
    
----------

lambda_1 : float
    Arrival rate of new critical patients (patients per day), mean is lambda_1 patients per day. Number of arrivals is a Poisson Process.
lambda_2 : float
    Arrival rate of new semi-critical patients (patients per day), mean is lambda_2 patients per day. Number of arrivals is a Poisson Process.
mu1 : float
    Parameter to determine the critical state treatment time. The treatment time is log-normally distributed. This is the mean of the underlying normal distribution.
sigma1 : float
    Parameter to determine the critical state treatment time. The treatment time is log-normally distributed. This is the standard deviation of the underlying normal distribution.
mu2 : float
    Parameter to determine the semi-critical state treatment time. The treatment time is log-normally distributed. This is the mean of the underlying normal distribution.
sigma2 : float
    Parameter to determine the semi-critical state treatment time. The treatment time is log-normally distributed. This is the standard deviation of the underlying normal distribution.
c1 : int
    Number of beds in the ICU.
c2 : int
    Number of beds in the SDU.
Q : int
    Maximum queue length. To represent the situation with no queues, set Q to be 1 and theta to be 0, then combine the balking proportion and abandoning proportion.
theta : float
    Rate of patients abandoning the queue, mean is theta days. Waiting times before abandoning are exponentially distributed.
p : float
    Probability of a critical patient immediately transitioning to the semi-critical state.
p_cd : float
    Probability of dying after leaving critical state.
p_cd_offplace : float
    Increase probability of dying after leaving critical state if offplaced.
p_cc_R : float
    Probability of a critical patient returning as a critical patient.
p_cc_R_offplace : float
    Increase probability of returning as critical after leaving critical state if offplaced.
p_csc_R : float
    Probability of a critical patient returning as a semi-critical patient.
p_csc_R_offplace : float
    Increase probability of returning as semi-critical after leaving critical state if offplaced.
p_scd : float
    Probability of dying after leaving semi-critical state.
p_scd_offplaced : float
    Increased probability of dying after leaving semi-critical state if offplaced.
p_scc_R : float
    Probability of a semi-critical patient returning as a critical patient.
p_scc_R_offplace : float
    Increase probability of returning as critical after leaving semi-critical state if offplaced.
p_scsc_R : float
    Probability of a semi-critical patient returning as a semi-critical patient.
p_scsc_R_offplace : float
    Increase probability of returning as semi-critical after leaving semi-critical state if offplaced.
    
delta : float
    Average waiting time before a returning patient returns. Waiting times are exponentially distributed.
x : float
    How much more time critical patients need if offplaced.
y : float
    How much more time semi-critical patients need if offplaced.
num_of_days : int
    Number of days to simulate. First 1/10 of the simulation is used for warm-up.
sampling_freq : float
    Sampling frequency to record the number of offplaced critical and semi-critical patients.

Returns\

-------

total_patients : int
    Total Number of Nurses involved.
total_critical_patients : int
    Total Number of Nurses that have ever been in the ciritical state (arrived as critical or returned as critical)
total_abandon_time : int
    Total number of times where a critical patient in a queue has abandoned waiting in the queue.
total_balked_time : int
    Total number of times where a new critical patient has been balked.
total_bumping_time : int
    Total number of times where a semi-critical patient has been bumped to the ward.
total_death_critical : int
    Total number of deaths after being in the critical state.
total_death_semi_critical : int
    Total number of deaths after being in the semi-critcal state.
total_offplacement_times_critical : float
    Total time (in days) that critical patients have been offplaced.
average_offplaced_critical : float
    Average number offplaced critcal patients at each timepoint sampled using sampling_freq.
total_offplacement_times_semi_critical : float
    Total time (in days) that semi-critical patients have been offplaced.
average_offplaced_semi_critical : float
    Average number offplaced semi-critcal patients at each timepoint sampled using sampling_freq.
total_queue_times : float
    Total time (in days) that critical patients have been in a queue.
queue_times_count : int
    Total number of times when a critical patient is in a queue or could have been in a queue.
"""


    event_queue = []
    heapq.heappush(event_queue, (0, "critical_arrival", 0))  # First critical arrival at time 0
    heapq.heappush(event_queue,(0.1, "semi_critical_arrival", 1))  # First semi-critical arrival at time 0.1
    # Tracking variables
    critical_customers = set()
    customers = {}
    total_patients = 0
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
    total_return_critical_critical = 0
    total_return_critical_semi_critical = 0
    total_return_semi_critical_critical = 0
    total_return_semi_critical_semi_critical = 0
    total_critical_treatments = 0
    total_semi_critical_treatments = 0

    event_time = 0
    recording_flag = False
    event_num = 0
    record_time_1 = num_of_days / 10
    record_time_2 = num_of_days / 10
    while event_time <= num_of_days:
        event_num += 1
        if event_num > 1000000:
            break
        
        if event_time >= num_of_days / 10:
            recording_flag = True

        # Process next event
        event_time, event_type, customer_id = heapq.heappop(event_queue)
        
       # Record number of patients in the ICU, SDU and the general ward
        while event_time >= record_time_1:
            if sampling_freq <= 0:
                break
            offplaced_occupancy_2_freq_1.append(sum(1 for server_status in servers_2 if server_status.get("kind") == "1"))
            offplaced_occupancy_3_freq_1.append(len(servers_3))
            time_points_freq_1.append(record_time_1)
            record_time_1 += sampling_freq
        
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



        if event_type == "critical_arrival":  # New critical patient
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

            if len(critical_queue) + len(return_critical_queue) >= Q:  # Queue length exceeded, balk
                customer_data["balked"] = True
                customers[customer_id] = customer_data
                if recording_flag:
                    total_balked_times += 1
                    total_patients += 1
                    critical_customers.add(customer_id)
                next_arrival_time = event_time + np.random.exponential(1 / lambda_1)
                heapq.heappush(event_queue, (next_arrival_time, "critical_arrival", customer_id + 2))
                continue
            else:
                critical_queue.append(customer_id)
                customers[customer_id] = customer_data
                if recording_flag:
                    total_patients += 1
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
            elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) < c2):  # No free bed in ICU but there are people to bump to the SDU
                service_start_time = event_time
                offplaced = False
                for server_status in servers_1[:]:
                    if server_status["kind"] == "2":
                        bumped_customer_id = server_status["customer_id"]
                        heapq.heappush(event_queue,(event_time, "bump_to_SDU", bumped_customer_id))
                        servers_1.remove(server_status)
                        break
            elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) >= c2):  # No free bed in ICU or SDU but there are people to bump to the ward
                service_start_time = event_time
                offplaced = False
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
                service_time1 = np.random.lognormal(mu1, sigma1)
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
                # print(f"updated patient {served_customer_id} information")

            else:  # Patient is offplaced
                queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
                service_time1 = np.random.lognormal(mu1, sigma1) * x
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


            # Schedule next arrival
            next_arrival_time = event_time + np.random.exponential(1 / lambda_1)
            heapq.heappush(event_queue, (next_arrival_time, "critical_arrival", customer_id + 2))



        
        elif event_type == "return_critical":  # Patient returning as critical
            customers[customer_id]["arrival_times"].append(event_time)
            return_critical_queue.append(customer_id)
            abandonment_time = event_time + np.random.exponential(1 / theta)
            critical_customers.add(customer_id)
            heapq.heappush(event_queue, (abandonment_time, "abandon", customer_id))
            served_customer_id = return_critical_queue[0]
            return_critical_queue.remove(served_customer_id)


            if len(servers_1) < c1:  # Free bed in ICU
                service_start_time = event_time
                offplaced = False
            elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) < c2):  # No free bed in ICU but there are people to bump to the SDU
                service_start_time = event_time
                offplaced = False
                for server_status in servers_1[:]:
                    if server_status["kind"] == "2":
                        bumped_customer_id = server_status["customer_id"]
                        heapq.heappush(event_queue,(event_time, "bump_to_SDU", bumped_customer_id))
                        servers_1.remove(server_status)
                        break
            elif (sum(1 for server_status in servers_1 if (server_status["kind"] == "1")) < c1) and (len(servers_2) >= c2):  # No free bed in ICU or SDU but there are people to bump to the ward
                service_start_time = event_time
                offplaced = False
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
            else:  # No way to begin treatment, put back to queue
                return_critical_queue.insert(0, served_customer_id)
                continue

            if not offplaced :  # Patient is not balked nor offplaced
                queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
                service_time1 = np.random.lognormal(mu1, sigma1)
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

            else:  # Patient is offplaced
                queue_time = event_time - customers[served_customer_id]["arrival_times"][-1]
                service_time1 = np.random.lognormal(mu1, sigma1) * y
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



        
        elif event_type == "semi_critical_arrival":  # New semi-critical patient
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
                total_patients += 1

            if len(servers_2) < c2:  # Free bed in SDU
                service_start_time = event_time
                offplaced = False
                service_time2 = np.random.lognormal(mu2, sigma2)
                customers[served_customer_id]["stage2_start_times"].append(service_start_time)
                customers[served_customer_id]["stage2_service_times"].append(service_time2)
                customers[served_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time2
                heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
                servers_2.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            elif len(servers_1) < c1:  # No Free bed in SDU, free bed in ICU
                service_start_time = event_time
                offplaced = False
                service_time2 = np.random.lognormal(mu2, sigma2)
                customers[served_customer_id]["stage2_start_times"].append(service_start_time)
                customers[served_customer_id]["stage2_service_times"].append(service_time2)
                customers[served_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time2
                heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
                servers_1.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            elif len(servers_1) >= c1 and len(servers_2) >= c2:  # No free bed in ICU and SDU, offplace
                service_start_time = event_time
                offplaced = True
                service_time2 = np.random.lognormal(mu2, sigma2)
                customers[served_customer_id]["stage2_start_times"].append(service_start_time)
                customers[served_customer_id]["stage2_service_times"].append(service_time2)
                customers[served_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time2
                if recording_flag:
                    total_offplacement_times_semi_critical += service_time2
                    offplacement_times_count_semi_critical += 1
                heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
                servers_3.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            # Schedule next arrival
            next_arrival_time = event_time + np.random.exponential(1 / lambda_2)
            heapq.heappush(event_queue, (next_arrival_time, "semi_critical_arrival", customer_id + 2))



        
        elif event_type == "return_semi_critical":  # Patient returning as semi-critical
            customers[customer_id]["arrival_times"].append(event_time)
            served_customer_id = customer_id

            if len(servers_2) < c2:  # Free bed in SDU
                service_start_time = event_time
                offplaced = False
                service_time2 = np.random.lognormal(mu2, sigma2)
                customers[served_customer_id]["stage2_start_times"].append(service_start_time)
                customers[served_customer_id]["stage2_service_times"].append(service_time2)
                customers[served_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time2
                heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
                servers_2.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            elif len(servers_1) < c1:  # No free bed in SDU, free bed in ICU
                service_start_time = event_time
                offplaced = False
                service_time2 = np.random.lognormal(mu2, sigma2)
                customers[served_customer_id]["stage2_start_times"].append(service_start_time)
                customers[served_customer_id]["stage2_service_times"].append(service_time2)
                customers[served_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time2
                heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
                servers_1.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})
            elif len(servers_1) >= c1 and len(servers_2) >= c2:  # No free bed in ICU and SDU, offplace
                service_start_time = event_time
                offplaced = True
                service_time2 = np.random.lognormal(mu2, sigma2) * y
                customers[served_customer_id]["stage2_start_times"].append(service_start_time)
                customers[served_customer_id]["stage2_service_times"].append(service_time2)
                customers[served_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time2
                if recording_flag:
                    total_offplacement_times_semi_critical += service_time2
                    offplacement_times_count_semi_critical += 1
                heapq.heappush(event_queue, (leave_time, "end_semi_critical", served_customer_id))
                servers_3.append({"time": leave_time, "kind": "2", "customer_id": served_customer_id})




        elif event_type == "critical_to_semi_critical":
            service_time2 = np.random.lognormal(mu2, sigma2)
            leave_time = service_time2 + event_time
            offplaced = False
            customers[customer_id]["stage2_start_times"].append(event_time)
            customers[customer_id]["stage2_service_times"].append(service_time2)
            customers[customer_id]["offplaced"].append(offplaced)
            for server_status in servers_1[:]:
                    if (server_status["customer_id"] == customer_id) and (server_status["kind"] == "1"):
                        servers_1.remove(server_status)
                        servers_1.append({"time": leave_time, "kind": "2", "customer_id": customer_id})
                        break
            for server_status in servers_2[:]:
                    if (server_status["customer_id"] == customer_id) and (server_status["kind"] == "1"):
                        servers_2.remove(server_status)
                        servers_2.append({"time": leave_time, "kind": "2", "customer_id": customer_id})
                        break
            heapq.heappush(event_queue, (leave_time, "end_semi_critical", customer_id))



        
        elif event_type == "abandon":  # Patient in queue abandon the queue
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



        
        elif event_type == "bump_to_SDU":  # Semi-critical patient is bumped to the SDU
            leave_time = customers[customer_id]["stage2_start_times"][-1] + customers[customer_id]["stage2_service_times"][-1]
            servers_2.append({"time": leave_time, "kind": 2, "customer_id": customer_id})



        
        elif event_type == "bump_to_ward":  # Semi-critical patient is bumped to the ward
            remaining_time = customers[customer_id]["stage2_start_times"][-1] + customers[customer_id]["stage2_service_times"][-1] - event_time
            leave_time = event_time + remaining_time * y
            servers_3.append({"time": leave_time, "kind": 2, "customer_id": customer_id})
            customers[customer_id]["bumped_times"] += 1
            customers[customer_id]["stage2_service_times"][-1] += remaining_time * (y - 1)
            if recording_flag:
                total_offplacement_times_semi_critical += remaining_time * y
                offplacement_times_count_semi_critical += 1
            heapq.heappush(event_queue, (event_time + remaining_time * y, "end_semi_critical", customer_id))



        
        elif event_type == "leave_ICU":  # Someone in the ICU finished treatment
            customer_found = False
            for server_status in servers_1[:]:
                if server_status["customer_id"] == customer_id:
                    servers_1.remove(server_status)
                    customer_found = True
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
                customers[transfer_customer_id]["stage1_service_times"][-1] = event_time + longest_remaining_time / x - customers[transfer_customer_id]["stage1_start_times"][-1]
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
                service_time1 = np.random.lognormal(mu1, sigma1)
                offplaced = False
                customers[new_customer_id]["stage1_start_times"].append(service_start_time)
                customers[new_customer_id]["queue_times"].append(queue_time)
                customers[new_customer_id]["stage1_service_times"].append(service_time1)
                customers[new_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time1
                heapq.heappush(event_queue, (leave_time, "end_critical", new_customer_id))
                servers_1.append({"time": leave_time, "kind": "1", "customer_id": new_customer_id})
            elif critical_queue:  # New critical patient in the queue
                new_customer_id = critical_queue[0]
                critical_queue.remove(new_customer_id)
                service_start_time = event_time
                queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
                if recording_flag:
                    total_queue_times += queue_time
                    queue_times_count += 1
                service_time1 = np.random.lognormal(mu1, sigma1)
                offplaced = False
                customers[new_customer_id]["stage1_start_times"].append(service_start_time)
                customers[new_customer_id]["queue_times"].append(queue_time)
                customers[new_customer_id]["stage1_service_times"].append(service_time1)
                customers[new_customer_id]["offplaced"].append(offplaced)
                leave_time = service_start_time + service_time1
                heapq.heappush(event_queue, (leave_time, "end_critical", new_customer_id))
                servers_1.append({"time": leave_time, "kind": "1", "customer_id": new_customer_id})
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



        
        elif event_type == "leave_SDU":  # Someone in the SDU finished treatment
            customer_found = False
            for server_status in servers_2[:]:
                if server_status["customer_id"] == customer_id:
                    servers_2.remove(server_status)
                    customer_found = True
                    break
            if not customer_found:
                continue
            if servers_3:  # Transfer offplaced semi_critical patient to SDU
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
                servers_2.append({"time": event_time + longest_remaining_time / y, "kind": "2", "customer_id": transfer_customer_id})
                customers[transfer_customer_id]["stage2_service_times"][-1] = event_time + longest_remaining_time / y - customers[transfer_customer_id]["stage2_start_times"][-1]
                heapq.heappush(event_queue, (event_time + longest_remaining_time / y, "end_semi_critical", transfer_customer_id))
            elif return_critical_queue:  # Offplace returned critical patients
                new_customer_id = return_critical_queue[0]
                return_critical_queue.remove(new_customer_id)
                service_start_time = event_time
                queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
                if recording_flag:
                    total_queue_times += queue_time
                    queue_times_count += 1
                service_time1 = np.random.lognormal(mu1, sigma1) * x
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
            elif critical_queue:  # Offplace new critical patients
                new_customer_id = critical_queue[0]
                critical_queue.remove(new_customer_id)
                service_start_time = event_time
                queue_time = service_start_time - customers[new_customer_id]["arrival_times"][-1]
                if recording_flag:
                    total_queue_times += queue_time
                    queue_times_count += 1
                service_time1 = np.random.lognormal(mu1, sigma1) * x
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


        
        elif event_type == "leave_ward":  # Someone in the ward finished treatment
            customer_found = False
            for server_status in servers_3[:]:
                if server_status["customer_id"] == customer_id:
                    servers_3.remove(server_status)
                    customer_found = True
                    if recording_flag:
                        total_bumping_times += 1
                    break



        
        elif event_type == "end_critical": # Critical patient finished critical treatment
            end_time_1 = 0
            end_time_2 = 0
            try:
                end_time_1 = customers[customer_id]["stage1_start_times"][-1] + customers[customer_id]["stage1_service_times"][-1]

            except Exception as e:
                pass

            if (abs(end_time_1 - event_time) < 1e-10) or (abs(end_time_2 - event_time) < 1e-10):
                if recording_flag:
                    total_critical_treatments += 1
                decision = np.random.rand()
                offplace_penalty = customers[customer_id]["offplaced"][-1]
                if decision < p:
                    # Schedule transition to semi-critical
                    add_returning_event_flag = True
                    for event in event_queue:
                        if (event[2] == customer_id) and ((event[1] == "return_critical") or (event[1] == "return_semi_critical") or (event[1] == "critical_to_semi_critical")):
                            add_returning_event_flag = False
                    if add_returning_event_flag:
                        heapq.heappush(event_queue, (event_time, "critical_to_semi_critical", customer_id))
                elif decision < p + p_cc_R + p_cc_R_offplace * offplace_penalty:
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
                        if recording_flag:
                            total_return_critical_critical += 1
                elif decision < p + p_cc_R + p_cc_R_offplace * offplace_penalty + p_csc_R + p_csc_R_offplace * offplace_penalty:
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
                        if recording_flag:
                            total_return_critical_semi_critical += 1
                elif decision < p + p_cc_R + p_cc_R_offplace * offplace_penalty + p_csc_R + p_csc_R_offplace * offplace_penalty + p_cd + p_cd_offplace * offplace_penalty:
                    customers[customer_id]["died"] = True
                    if recording_flag:
                        total_death_critical += 1
                    heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                    heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))
                else:
                    # Schedule leave without returning to the system
                    heapq.heappush(event_queue,(event_time, "leave_SDU", customer_id))
                    heapq.heappush(event_queue,(event_time, "leave_ICU", customer_id))



        
        elif event_type == "end_semi_critical":  # Semi-critical patient finished semi-critical treatment
            end_time_1 = 0
            end_time_2 = 0
            try:
                end_time_2 = customers[customer_id]["stage2_start_times"][-1] + customers[customer_id]["stage2_service_times"][-1]

            except Exception as e:
                pass

            if (abs(end_time_1 - event_time) < 1e-10) or (abs(end_time_2 - event_time) < 1e-10):
                if recording_flag:
                    total_semi_critical_treatments += 1
                decision = np.random.rand()
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
                        if recording_flag:
                            total_return_semi_critical_semi_critical += 1
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
                        if recording_flag:
                            total_return_semi_critical_critical += 1
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

    total_critical_patients = len(critical_customers)
    average_offplaced_critical = np.mean(offplaced_occupancy_2_freq_1)
    average_offplaced_semi_critical = np.mean(offplaced_occupancy_3_freq_1)

 
    return (total_patients, #0
            total_critical_patients, #1
            total_abandon_times, #2
            total_balked_times, #3
            total_bumping_times, #4
            total_death_critical, #5
            total_death_semi_critical, #6 
            total_return_critical_critical, #7
            total_return_critical_semi_critical, #8
            total_return_semi_critical_critical, #9
            total_return_semi_critical_semi_critical, #10 
            total_offplacement_times_critical, #11
            average_offplaced_critical, #12
            total_offplacement_times_semi_critical, #13
            average_offplaced_semi_critical, #14
            total_critical_treatments, #15
            total_semi_critical_treatments, #16
            total_queue_times, #17
            queue_times_count) #18
