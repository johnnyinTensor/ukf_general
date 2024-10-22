#include "ukf_10m.h"


int main(){
    // ... (initialization code) ...
    // Get sensor data (implementation depends on your hardware)
    float sx, sy, sz, gx, gy, gz, mx, my, mz, sx_r, sy_r, sz_r, mx_r, my_r, mz_r;
    float dt = 0.1;  // Get time since last iteration
    double moi[3][3] = {{0.0333, 0, 0}, 
                          {0, 0.0060, 0}, 
                          {0, 0, 0.0333}};
    ukf_init((const double*)moi);

    while (1) {

        // Clear previous sensor data
        ukf_sensor_clear();

        // Set new sensor data
        ukf_sensor_set_fss(sx, sy, sz);
        ukf_sensor_set_gyroscope(gx, gy, gz);
        ukf_sensor_set_magnetometer(mx, my, mz);

        // Set reference data
        ukf_set_sun_ref(sx_r, sy_r, sz_r);
        ukf_set_mag_ref(mx_r, my_r, mz_r);

        // Iterate the filter
        ukf_iterate(dt);

        // Get the updated state
        struct ukf_state_t state;
        ukf_get_state(&state);
    }

    return 0;
}
