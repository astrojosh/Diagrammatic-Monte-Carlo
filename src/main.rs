use rand::distributions::{Distribution, Uniform};
use std::time::{Instant};

const U_0: f64 = 1.0;
const R_0: f64 = 1.0;
const Q_0: f64 = 1.0;
const Q_MAX: f64 = 100.0;
const N_BINS: usize = 1000;
const D_SIGMA: f64 = 2.0;

// static DELTA_Q = Q_MAX / (N_BINS as f64);
static DELTA_Q: f64 = 0.1;

const EPSILON: f64 = 1e-7;
const LOWEST_ORDER: u32 = 0;
const HIGHEST_ORDER: u32 = 2;
const PI: f64 = 3.14159265358979323846;

fn main() {

    // Start timer to measure run time
    let start = Instant::now();

    let answer = 0.371817;

    let mut step = 0;
    let mut s = 1.0;
	let mut rejects = 0;
	let mut hist = [0.0; N_BINS].to_vec();
	let mut zsigma = 0;
	let mut order = 0;
    let mut chi = 0.0;
    let mut qp = 0.0;
    let mut q = 0.0;
    let uniform_range = Uniform::from(0.0..1.0);
    let mut rng = rand::thread_rng();

    while _ in 0..N {
        step += 1;

        if order == 0 {
            zsigma += 1;
        }
        else if q < Q_MAX {
            let index = (q / DELTA_Q) as usize;
            hist[index] += s;
        }

        if LOWEST_ORDER == order {

            let srcOrder = order;
            let destOrder = order + 1;

            let mut news = 0.0;
            let r = uniform_range.sample(&mut rng);
            let mut ratio = 1.0;
            let high_order = srcOrder.max(destOrder);
            let low_order = srcOrder.min(destOrder);
        
            if destOrder > srcOrder {

                let r = uniform_range.sample(&mut rng);
            
                if destOrder == 1 {
                    q = Q_0 * r / (Q_0 / Q_MAX + 1.0 - r);
                }
            }
        
            let higher = Level(high_order, q, qp, chi, &mut hist, zsigma);
            let lower = Level(low_order, q, qp, chi, &mut hist, zsigma);
        
            if destOrder < srcOrder {
                ratio = (lower / higher).abs();
                news = lower.signum();
            }
            else {
                ratio = (higher / lower).abs();
                news = higher.signum();
            }
        
            if ratio < 1.0 && r > ratio {
                rejects += 1;
            }
            else {
                s = news;
                order = destOrder;
            }
        }
        else if HIGHEST_ORDER == order {
            // Update(order, order - 1);

            let srcOrder = order;
            let destOrder = order - 1;

            let mut news = 0.0;
            let r = uniform_range.sample(&mut rng);
            let mut ratio = 1.0;
            let high_order = srcOrder.max(destOrder);
            let low_order = srcOrder.min(destOrder);

            let higher = Level(high_order, q, qp, chi, &mut hist, zsigma);
            let lower = Level(low_order, q, qp, chi, &mut hist, zsigma);

            if destOrder < srcOrder {
                ratio = (lower / higher).abs();
                news = lower.signum();
            }

            if ratio < 1.0 && r > ratio {
                rejects += 1;
            }
            else {
                s = news;
                order = destOrder;
            }
        }
        else {
            let r = uniform_range.sample(&mut rng);
            if r > 0.5 {
                // Update(order, order + 1);
                let srcOrder = order;
                let destOrder = order + 1;
    
                let mut news = 0.0;
                let r = uniform_range.sample(&mut rng);
                let mut ratio = 1.0;
                let high_order = srcOrder.max(destOrder);
                let low_order = srcOrder.min(destOrder);
            
                if destOrder > srcOrder {
                    // Variables(destOrder);
    
                    let uniform_range = Uniform::from(0.0..1.0);
                    let r = uniform_range.sample(&mut rng);
                
                    if destOrder == 1 {
                        q = Q_0 * r / (Q_0 / Q_MAX + 1.0 - r);
                    }
                    else if destOrder == 2 {
                        qp = Q_0 * r / (1.0 - r);
                        let r = uniform_range.sample(&mut rng);
                        chi = 2.0 * r - 1.0;
                    }
                }
            
                let higher = Level(high_order, q, qp, chi, &mut hist, zsigma);
                let lower = Level(low_order, q, qp, chi, &mut hist, zsigma);
            
                if destOrder < srcOrder {
                    ratio = (lower / higher).abs();
                    news = lower.signum();
                }
                else {
                    ratio = (higher / lower).abs();
                    news = higher.signum();
                }
            
                if ratio < 1.0 && r > ratio {
                    rejects += 1;
                }
                else {
                    s = news;
                    order = destOrder;
                }
            }
            else {
                // Update(order, order - 1);
                let srcOrder = order;
                let destOrder = order - 1;
    
                let mut news = 0.0;
                let r = uniform_range.sample(&mut rng);
                let mut ratio = 1.0;
                let high_order = srcOrder.max(destOrder);
                let low_order = srcOrder.min(destOrder);
    
                let higher = Level(high_order, q, qp, chi, &mut hist, zsigma);
                let lower = Level(low_order, q, qp, chi, &mut hist, zsigma);

                if destOrder < srcOrder {
                    ratio = (lower / higher).abs();
                    news = lower.signum();
                }

                if ratio < 1.0 && r > ratio {
                    rejects += 1;
                }
                else {
                    s = news;
                    order = destOrder;
                }
            }
        }

        if step % 10000000 == 0 {
            // Print();
            println!("Rejected diagram changes  =  {}\n", rejects);
            println!("Normalization statistics  =  {}\n", zsigma);
            println!("MC -f(0) value            =  {}\n", -f(0.0, &mut hist, zsigma));
            println!("Exact answer              =  {}\n", answer);
            println!("-------------------------------------\n");
        }
    }

    // ----------------------------------------------

    // End timer to measure run time
    let duration = start.elapsed();

    println!("\n--------------------------------------------------------------------");
    // Print run times
    println!("Total time\t\t\t= {:?}", duration);
    println!("--------------------------------------------------------------------\n");
}

fn U(q: f64) -> f64 {
    if q < EPSILON {
        2.0 * U_0 * f64::powi(R_0, 3) / 3.0
    }
    else {
        2.0 * U_0 * ((q * R_0).sin() - q * R_0 * (q * R_0).cos()) / f64::powi(q, 3)
    }
}

fn f(q: f64, hist: &mut Vec<f64>, zsigma: u32) -> f64 {
    if q > Q_MAX || zsigma < 100000 {
        -U(q)
    }
    else {
        let index = (q / DELTA_Q) as usize;
        hist[index] * D_SIGMA / ((zsigma as f64) * DELTA_Q)
    }
}

fn WQp(q: f64) -> f64 {
    Q_0 / f64::powi(Q_0 + q, 2)
}

fn WQ(q: f64) -> f64 {
    ((Q_0 + Q_MAX) / Q_MAX) * WQp(q)
}

fn WChi(chi: f64) -> f64 {
    0.5
}

fn Level(order: u32, q: f64, qp: f64, chi: f64, hist: &mut Vec<f64>, zsigma: u32) -> f64 {
    let mut value = 0.0;
    let mut norm = 1.0;

    if order == 0 {
        value = D_SIGMA;
    }
    else if order == 1 {
        norm = 0.5;
        value = -U(q) / WQ(q)
    }
    else if order == 2 {
        value = (f64::powi(q, 2) + f64::powi(qp, 2) - 2.0 * q * qp * chi).sqrt();
        value = -U(value) * f(qp, hist, zsigma) / (PI * WQp(qp) * WChi(chi) * WQ(q));
    }

    norm * value
}