use lazy_static::lazy_static;
use reqwest_middleware::{ClientBuilder, ClientWithMiddleware};
use reqwest_retry::{policies::ExponentialBackoff, RetryTransientMiddleware};
use tokio::time::Duration;

lazy_static! {
    static ref CLIENT: ClientWithMiddleware = {
        let retry_policy = ExponentialBackoff::builder().build_with_max_retries(3);
        let client = reqwest::Client::builder()
            // .timeout(Duration::from_secs(300))
            .connect_timeout(Duration::from_secs(20))
            // .tcp_keepalive(Duration::from_secs(300))
            .build()
            .expect("reqwest::Client::new()");

        ClientBuilder::new(client)
            .with(RetryTransientMiddleware::new_with_policy(retry_policy))
            .build()
    };
}

pub fn retry_client() -> &'static ClientWithMiddleware {
    &CLIENT
}
